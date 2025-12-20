import arxiv
from datetime import datetime, timedelta, timezone

def search_arxiv(query, max_results=30):
    """
    Searches ArXiv for papers matching the query from the last 7 days.
    """
    client = arxiv.Client()
    
    # Construct search query (Sorted by Submitted Date)
    search = arxiv.Search(
        query=query,
        max_results=max_results,
        sort_by=arxiv.SortCriterion.SubmittedDate
    )
    
    # Calculate cutoff date (7 days ago) in UTC
    seven_days_ago = datetime.now(timezone.utc) - timedelta(days=7)

    papers = []
    
    try:
        # Iterate through search results
        for result in client.results(search):
            
            # Optimization: Since results are sorted by date (newest first),
            # once we hit a paper older than 7 days, we can stop the loop entirely.
            if result.published < seven_days_ago:
                break 

            journal_name = result.journal_ref if result.journal_ref else "ArXiv Preprint"
            
            # Format data into a dictionary
            papers.append({
                "source": "ArXiv",
                "journal": journal_name,
                "title": result.title,
                "url": result.entry_id,
                "abstract": result.summary.replace("\n", " "), # Remove newlines for cleaner summary
                "date": result.published.strftime('%Y-%m-%d')
            })
            
        return papers

    except Exception as e:
        print(f"❌ ArXiv Search Error: {e}")
        return []