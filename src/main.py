import os
import sys
import time
import argparse
import json

# Add current directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from dotenv import load_dotenv
from src.search_arxiv import search_arxiv
from src.search_pubmed import search_pubmed
from src.summarizer import summarize_abstract
from src.notifier import send_to_discord

load_dotenv()

# === Default Configuration (Fallback) ===
DEFAULT_CONFIG = [
    {
        "source": "pubmed",
        "query": '"Chromatin" OR "Hi-C"',
        "webhook_env": "DISCORD_WEBHOOK_URL"
    }
]

# === Webhook Shortcut Map ===
# Maps short codes in the text file to actual environment variable names
WEBHOOK_MAP = {
    "DEFAULT": "DISCORD_WEBHOOK_URL",
    # Optional
    "CUSTOM1": "DISCORD_WEBHOOK_URL_CUSTOM_1",
    "CUSTOM2": "DISCORD_WEBHOOK_URL_CUSTOM_2"
}

def load_config_from_file(filepath):
    """Loads search configuration from a text file."""
    if not os.path.exists(filepath):
        print(f"ℹ️ Config file not found at: {filepath}")
        print("   -> Using internal DEFAULT configuration.")
        return DEFAULT_CONFIG

    print(f"📂 Loading config from: {filepath}")
    config_list = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line_no, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = [p.strip() for p in line.split('|')]
                
                if len(parts) < 2:
                    print(f"⚠️ Warning (Line {line_no}): Invalid format. Skipping.")
                    continue
                
                source = parts[0].lower()
                query = parts[1]
                webhook_key = "DEFAULT"
                
                if len(parts) >= 3 and parts[2]:
                    webhook_key = parts[2].upper()
                
                env_var_name = WEBHOOK_MAP.get(webhook_key, webhook_key)
                
                config_list.append({
                    "source": source,
                    "query": query,
                    "webhook_env": env_var_name
                })
                
        if not config_list:
            print("⚠️ Config file is empty. Using default.")
            return DEFAULT_CONFIG
            
        print(f"✅ Loaded {len(config_list)} search rules.")
        return config_list

    except Exception as e:
        print(f"❌ Error reading config file: {e}")
        return DEFAULT_CONFIG

def main():
    parser = argparse.ArgumentParser(description="Paper Digest: Automated Paper Collector & Summarizer")
    parser.add_argument("--config", default="search_config.txt", help="Path to configuration file")

    parser.add_argument("--test", action="store_true", help="Run in test mode (process only 1 paper per rule)")
    args = parser.parse_args()

    print("=== Paper Digest Started ===")
    if args.test:
        print("🧪 TEST MODE ACTIVE: Processing only 1 paper per search rule.")
    
    search_config = load_config_from_file(args.config)
    
    for config in search_config:
        print(f"\n🔍 Searching {config['source']} for: {config['query']}")
        
        webhook_url = os.getenv(config['webhook_env'])
        
        if not webhook_url:
            print(f"⚠️ Webhook '{config['webhook_env']}' not found. Trying default webhook.")
            webhook_url = os.getenv("DISCORD_WEBHOOK_URL")
        
        if not webhook_url:
             print("❌ Error: No webhook URL found. Skipping this rule.")
             continue

        papers = []
        if config['source'] == "arxiv":
            papers = search_arxiv(config['query'])
        elif config['source'] == "pubmed":
            papers = search_pubmed(config['query'])
            
        if not papers:
            print(f"No new papers found.")
            continue

        if args.test:
            print(f"🧪 Test mode: Limiting results to 1 out of {len(papers)} found.")
            papers = papers[:1]
            
        for paper in papers:
            print(f"Summarizing: {paper['title']}")
            
            summary = ""
            max_retries = 3
            
            for attempt in range(max_retries):
                try:
                    summary = summarize_abstract(paper['abstract'])
                    
                    if "429" in summary or "Quota exceeded" in summary:
                        raise Exception("API Rate Limit Hit (Response content)")
                    
                    break # Success
                    
                except Exception as e:
                    print(f"⚠️ API Error (Attempt {attempt+1}/{max_retries}): {e}")
                    if attempt < max_retries - 1:
                        wait_time = 60
                        print(f"Waiting {wait_time} seconds before retrying...")
                        time.sleep(wait_time)
                    else:
                        summary = "⚠️ API Error: Could not generate summary due to rate limits."

            send_to_discord(paper, summary, webhook_url)

            if not args.test and len(papers) > 1:

                print("Sleeping for 30 seconds...")
                time.sleep(30)
            
    print("\n=== Finished successfully ===")

if __name__ == "__main__":
    main()