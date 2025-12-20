import os
import requests

def send_to_discord(paper_data, summary, webhook_url):
    """
    指定されたDiscord Webhookへメッセージを送信する関数
    """
    if not webhook_url:
            print("Error: Webhook URL is missing.")
            return

    # メッセージの整形
    content = f"""
**【新着論文: {paper_data['source']}】**
**掲載誌:** {paper_data['journal']}
**発行日:** {paper_data['date']}
**タイトル:** {paper_data['title']}
**リンク:** {paper_data['url']}

**[AI要約]**
{summary}

-----------------------------------
"""
    
    # 送信処理
    try:
        response = requests.post(webhook_url, json={"content": content})
        if response.status_code == 204:
            print(f"Sent: {paper_data['title']}")
        else:
            print(f"Discord Error: {response.status_code}")
    except Exception as e:
        print(f"Discord Connection Error: {e}")