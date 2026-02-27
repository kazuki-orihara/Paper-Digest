from Bio import Entrez
import os

# Set email from environment variable to identify the user to NCBI
# (Required by NCBI terms of use)
Entrez.email = os.getenv("ENTREZ_EMAIL", "paper_digest_bot@example.com")

def search_pubmed(keyword, max_results=1):
    """
    Search PubMed for papers matching the keyword within the last 7 days.
    """
    try:
        # 1. Search for IDs
        # 'reldate=7' limits search to the last 7 days
        # 'datetype="pdat"' specifies Publication Date
        handle = Entrez.esearch(
            db="pubmed", 
            term=keyword, 
            retmax=max_results, 
            sort="date", 
            reldate=7, 
            datetype="pdat"
        )
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list:
            return []

        # 2. Fetch Detailed Info
        handle = Entrez.efetch(
            db="pubmed", 
            id=id_list, 
            rettype="medline", 
            retmode="xml"
        )
        records = Entrez.read(handle)
        
        papers = []
        for article in records['PubmedArticle']:
            medline = article['MedlineCitation']['Article']
            
            # Extract Abstract
            abstract_text = ""
            if 'Abstract' in medline and 'AbstractText' in medline['Abstract']:
                abstract_list = medline['Abstract']['AbstractText']
                # Join paragraphs if abstract is a list
                if isinstance(abstract_list, list):
                    # Handle cases where list items might be objects with attributes
                    abstract_text = " ".join([str(item) for item in abstract_list])
                else:
                    abstract_text = str(abstract_list)

            # Skip papers without abstracts
            if abstract_text:
                # --- Extract Publication Date ---
                pub_date_data = medline.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                date_str = "Unknown"
                
                if 'Year' in pub_date_data:
                    date_str = pub_date_data['Year']
                    if 'Month' in pub_date_data:
                        date_str += f"-{pub_date_data['Month']}"
                    if 'Day' in pub_date_data:
                        date_str += f"-{pub_date_data['Day']}"
                elif 'MedlineDate' in pub_date_data:
                    date_str = pub_date_data['MedlineDate'] # e.g., "2024 Jan-Feb"

                # --- Extract Journal Name ---
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
        print(f"❌ PubMed Search Error: {e}")
        return []
