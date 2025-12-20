import os
import sys
import time
# カレントディレクトリをパスに追加（モジュール読み込みのため）
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from dotenv import load_dotenv
from src.search_arxiv import search_arxiv
from src.search_pubmed import search_pubmed
from src.summarizer import summarize_abstract
from src.notifier import send_to_discord

# .envファイルを読み込み（ローカル実行時のみ有効）
load_dotenv()

# === 設定エリア ===
# ここに検索したいキーワードを設定します
SEARCH_CONFIG = [
    {"source": "pubmed", # クロマチン構造関連
     "query": '"Chromatin" OR ("Hi-C" OR "TADs" OR "3D Genome")',
     "webhook_env":"DISCORD_WEBHOOK_URL_INT"
     },
# 1. 【超トップジャーナル】絶対に逃したくないやつ
    # Nature, Science, Cell, Molecular Cell, Nature Genetics など
    # 掲載数が比較的少ないので、埋もれにくい
    {
        "source": "pubmed",
        "query": '("Nature"[Journal] OR "Science"[Journal] OR "Cell"[Journal] OR "Molecular cell"[Journal] OR "Nature genetics"[Journal] OR "Nature cell biology"[Journal] OR "Genome biology"[Journal]) AND ("Chromatin" OR "Hi-C" OR "TADs" OR "3D Genome" OR "Nucleosome" OR "Epigenetics" OR "Phase separation")',
        "webhook_env": "DISCORD_WEBHOOK_URL_TOP"
    },

    # 2. 【メガジャーナル・特定分野誌】数は多いが重要なやつ
    # Nature Communications, Bioinformatics, NARなど
    # 分けておくことで、ここが大量にヒットしても上記1番を邪魔しない
    {
        "source": "pubmed",
        "query": '("Nature communications"[Journal] OR "Bioinformatics"[Journal] OR "Nucleic acids research"[Journal]) AND ("Chromatin" OR "Hi-C" OR "TADs" OR "3D Genome" OR "Nucleosome" OR "Epigenetics" OR "Phase separation")',
        "webhook_env": "DISCORD_WEBHOOK_URL_TOP" # 同じチャンネルに通知してOK
    },
    {"source": "arxiv", # from arxiv
     "query": '(cat:q-bio.BM OR cat:q-bio.GN) AND ("Chromatin" OR "Hi-C" OR "TADs" OR "3D Genome" OR "Nucleosome" OR "Phase separation" OR "Condensates")',
     "webhook_env":"DISCORD_WEBHOOK_URL_ARXIV"
     },
# 1. 【クロマチン物理】紐としての性質・アクティブマター
    # クロマチンの「動き」や「硬さ」を物理的に解析したもの
    {
        "source": "pubmed",
        "query": '("Chromatin"[Title/Abstract]) AND ("Polymer"[Title/Abstract] OR "Active matter"[Title/Abstract] OR "Viscoelasticity"[Title/Abstract] OR "Persistence length"[Title/Abstract] OR "Radius of gyration"[Title/Abstract])',
        "webhook_env": "DISCORD_WEBHOOK_URL_BIOPHYSICS"
    },

    # 2. 【相分離の物理】現象論ではなくメカニズム
    # 単に「相分離した」だけでなく、その物理的駆動力を議論しているもの
    {
        "source": "pubmed",
        "query": '("Phase separation"[Title/Abstract] OR "Condensates"[Title/Abstract]) AND ("Nucleation" OR "Spinodal" OR "Surface tension" OR "Ostwald ripening" OR "Material properties")',
        "webhook_env": "DISCORD_WEBHOOK_URL_BIOPHYSICS"
    },

    # 3. 【シミュレーション・計算】
    # クロマチンやRNAの挙動を計算機で再現したもの
    {
        "source": "pubmed",
        "query": '("Chromatin" OR "RNA" OR "Protein") AND ("Molecular dynamics"[Title/Abstract] OR "Coarse-grained"[Title/Abstract] OR "Monte Carlo simulation"[Title/Abstract] OR "Free energy landscape"[Title/Abstract])',
        "webhook_env": "DISCORD_WEBHOOK_URL_BIOPHYSICS"
    },
    
    # 4. 【一分子計測・力学】
    # 光ピンセットなどで実際に引っ張ったり力を測ったりしたもの
    {
        "source": "pubmed",
        "query": '("Chromatin" OR "DNA" OR "RNA") AND ("Optical tweezers" OR "Magnetic tweezers" OR "Force spectroscopy" OR "Stiffness")',
        "webhook_env": "DISCORD_WEBHOOK_URL_BIOPHYSICS"
    }
]
# =================

def main():
    print("=== Paper Notifier Started ===")
    
    for config in SEARCH_CONFIG:
        print(f"Searching {config['source']} for: {config['query']}")
        # 【重要修正】ここで環境変数名から「実際のURL」を取り出す！
        webhook_url = os.getenv(config['webhook_env'])
        
        # URLが設定されていない（.envに書き忘れている）場合はスキップしてエラーを防ぐ
        if not webhook_url:
            print(f"⚠️ Skipping: Webhook URL not found for variable '{config['webhook_env']}'")
            continue
        papers = []
        
        if config['source'] == "arxiv":
            papers = search_arxiv(config['query'])
        elif config['source'] == "pubmed":
            papers = search_pubmed(config['query'])
            
        if not papers:
            print(f"No papers found for {config['query']}")
            continue
            
        for paper in papers:
            print(f"Summarizing: {paper['title']}")
            summary = summarize_abstract(paper['abstract'])
            send_to_discord(paper, summary, webhook_url)

            print("Sleeping for 60 seconds to respect API limits...")
            time.sleep(60)  # API制限を考慮して60秒待機
            
    print("=== Use NotebookLM for deeper reading! ===")
    print("=== Finished ===")

if __name__ == "__main__":
    main()