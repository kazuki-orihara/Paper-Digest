# Paper Digest: Automated AI Paper Collector

Paper Digest is an automated tool that fetches the latest research papers from PubMed and ArXiv based on specific keywords, generates structured summaries using Google Gemini AI, and notifies your Discord server.

Automate your daily literature review and ensure you never miss important research.

## 🚀 Features

Automated Search: Supports both PubMed and ArXiv.

AI Summarization: Uses Gemini 2.5 Flash to generate professional summaries from abstracts.

Output is structured into 4 key sections: Background, Methods, Results, and Significance.

Discord Notifications: Sends clean, readable cards (Embeds) with the paper title, link, and AI summary.

Flexible Configuration: Easily modify search queries and destination channels by editing the search_config.txt file.

Smart Categorization: Route papers to different Discord channels (Webhooks) based on themes (e.g., "Biology", "Physics", "General").

## 🛠️ Setup Guide

### 1. Prerequisites

* `Discord Webhook URL`: Obtain this from the channel settings in your Discord server.

* `Google Gemini API Key`: Get it for free at Google AI Studio.

* `Python 3.10+`

### 2. Installation

Clone the repository and install the required libraries.

```sh
git clone https://github.com/citrusml/Paper-Digest.git
cd Paper-Digest
pip install -r requirements.txt
```

### 3. Environment Variables

Create a `.env` file in the root directory and set your API key and Webhooks.

```sh
# .env file example

# REQUIRED: Gemini API Key
GEMINI_API_KEY="your_google_api_key_here"

# REQUIRED: Default Discord Webhook URL
DISCORD_WEBHOOK_URL="https://discord.com/api/webhooks/xxxx/xxxx"

# OPTIONAL: If you want to separate channels (Specify BIO, PHYSICS, etc. in config file)
# DISCORD_WEBHOOK_URL_CUSTOM_1="your_webhook_url_for_custom1"
# DISCORD_WEBHOOK_URL_CUSTOM_2="your_webhook_url_for_custom2"

# Optional: Recommended for PubMed queries
ENTREZ_EMAIL="your_email@example.com"
```

### 4. Search Configuration (search_config.txt)

Edit `search_config.txt` to define the papers you want to collect. The format is very simple. (See search_config.txt for examples).

### 5. Usage

Run the following command to search, summarize, and notify:

```sh
python src/main.py
```

## 🤖 Automation (GitHub Actions)

This repository includes GitHub Actions configuration for automatic execution. Go to your repository's `Settings` > `Secrets and variables` > `Actions` and add the following secrets to run the bot automatically (default: every Monday).

`GEMINI_API_KEY`
`DISCORD_WEBHOOK_URL`