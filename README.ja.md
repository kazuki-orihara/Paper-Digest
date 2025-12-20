# Paper Digest: AI論文収集・要約ボット
Paper Digest は、PubMed や ArXiv から特定のキーワードで最新の論文を自動収集し、Google Gemini AI を用いて日本語で要約を作成、Discord に通知するツールです。

日々の論文チェックを自動化し、重要な研究を見逃さないようにサポートします。

## 🚀 できること (Features)
自動検索: PubMed と ArXiv の両方に対応。

AI要約: Gemini 1.5 Flash を使用し、アブストラクトを専門家向けの日本語に要約。

「背景」「手法」「結果」「意義」の4点に構造化して出力します。

Discord通知: 論文タイトル、リンク、要約を見やすいカード形式（Embed）で通知。

柔軟な設定: `search_config.txt` ファイルを編集するだけで、検索ワードや通知先を自由に変更可能。

カテゴリ分け: 「生物系」「物理系」「その他」など、テーマごとに通知するDiscordチャンネル（Webhook）を振り分け可能。

## 🛠️ 使い方の流れ (Setup Guide)

### 1. 準備するもの

* `Discord Webhook URL`: 通知を送りたいチャンネルの設定から取得します。

* `Google Gemini API Key`: Google AI Studio で無料で取得できます。

* `Python 3.10+`

### 2. インストール

リポジトリをクローンし、ライブラリをインストールします。

    ```sh
    git clone https://github.com/citrusml/Paper-Digest.git
    cd Paper-Digest
    pip install -r requirements.txt
    ```

### 3. 環境変数の設定:

プロジェクトのルートに `.env` ファイルを作成し、APIキーを設定します。


```sh
# .env ファイルの例

# 必須: Gemini APIキー
GEMINI_API_KEY="your_google_api_key_here"

# 必須: デフォルトのDiscord Webhook URL
DISCORD_WEBHOOK_URL="https://discord.com/api/webhooks/xxxx/xxxx"

# オプション: チャンネルを分けたい場合（設定ファイルで BIO, PHYSICS などと指定可能）
# DISCORD_WEBHOOK_URL_BIO="別のWebhookURL"
# DISCORD_WEBHOOK_URL_PHYSICS="別のWebhookURL"

# 推奨: PubMed queriesのため
ENTREZ_EMAIL="your_email@example.com"
```

### 4. 検索条件の設定 (search_config.txt)

`search_config.txt` を編集して、収集したい論文の条件を指定します。 書き方は非常にシンプルです。

### 5. 実行

以下のコマンドで検索・要約・通知を実行します。

    ```sh
    python src/main.py
    ```

## 🤖 自動化 (GitHub Actions)

このリポジトリには GitHub Actions の設定が含まれています。 リポジトリの `Settings` > `Secrets and variables` > `Actions` に以下の変数を登録すれば、決まった時間（デフォルトでは毎月曜日）に自動で実行されます。

GEMINI_API_KEY
DISCORD_WEBHOOK_URL