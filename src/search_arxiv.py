import arxiv
from datetime import datetime, timedelta, timezone

def search_arxiv(keyword, max_results=10):
    """
    ArXivから指定キーワードで論文を検索する関数
    """
    client = arxiv.Client()
    
    # 検索クエリの構築（提出順にソート）
    search = arxiv.Search(
        query=keyword,
        max_results=max_results,
        sort_by=arxiv.SortCriterion.SubmittedDate
    )
    
    # 7日前の日時を計算（UTCタイムゾーン）
    seven_days_ago = datetime.now(timezone.utc) - timedelta(days=7)

    papers = []
    for result in client.results(search):
            # 【判定】論文の発行日が、7日前より古い場合はスキップ（または終了）
            if result.published < seven_days_ago:
                # ArXivは新しい順に出てくるので、古いのが出たらそこで打ち切っても良いが
                #念のためcontinueにしておく
                continue

            journal_name = result.journal_ref if result.journal_ref else "ArXiv Preprint"
            # 必要情報を辞書にまとめる
            papers.append({
                "source": "ArXiv",
                "journal": journal_name,
                "title": result.title,
                "url": result.entry_id,
                "abstract": result.summary.replace("\n", " "), # 改行を除去
                "date": result.published.strftime('%Y-%m-%d')
            })
        
    return papers