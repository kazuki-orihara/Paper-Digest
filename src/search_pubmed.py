from Bio import Entrez
import datetime

# 自分のメールアドレスを設定（PubMedのルールとして必要）
Entrez.email = "yishuxiangyuan73@gmail.com" 

def search_pubmed(keyword, max_results=10):
    """PubMedからキーワードで論文を検索"""
    # 検索クエリの作成（日付指定などはここで調整可能）
    # ここではシンプルにキーワード検索して直近を取得
    
    try:
        # 1. IDを検索
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, sort="date",reldate=7,datetype="pdat")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list:
            return []

        # 2. 詳細情報を取得
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        
        papers = []
        for article in records['PubmedArticle']:
            medline = article['MedlineCitation']['Article']
            
            # アブストラクトがある場合のみ取得
            abstract_text = ""
            if 'Abstract' in medline and 'AbstractText' in medline['Abstract']:
                abstract_list = medline['Abstract']['AbstractText']
                # 複数の段落に分かれている場合があるため結合
                abstract_text = " ".join(abstract_list) if isinstance(abstract_list, list) else abstract_list

            if abstract_text:
                # --- 日付取得ロジックの追加 ---
                pub_date_data = medline.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                date_str = "Unknown"
                
                if 'Year' in pub_date_data:
                    date_str = pub_date_data['Year']
                    if 'Month' in pub_date_data:
                        date_str += f"-{pub_date_data['Month']}"
                    if 'Day' in pub_date_data:
                        date_str += f"-{pub_date_data['Day']}"
                elif 'MedlineDate' in pub_date_data:
                    date_str = pub_date_data['MedlineDate'] # "2024 Jan-Feb" などの形式

                # --- 雑誌名取得 ---
                journal_title = medline.get('Journal', {}).get('ISOAbbreviation', "")
                if not journal_title:
                    journal_title = medline.get('Journal', {}).get('Title', "Unknown Journal")

                pmid = article['MedlineCitation']['PMID']
                papers.append({
                    "source": "PubMed",
                    "journal": journal_title,
                    "title": medline['ArticleTitle'],
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "abstract": abstract_text,
                    "date": date_str
                })
        return papers
        
    except Exception as e:
        print(f"PubMed Error: {e}")
        return []