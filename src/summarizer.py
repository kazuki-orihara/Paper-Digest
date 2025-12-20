import os
import google.generativeai as genai

def summarize_abstract(text):
    """
    Geminiを使ってアブストラクトを要約する関数
    """
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        return "Error: API Key not found."

    genai.configure(api_key=api_key)
    model = genai.GenerativeModel('gemini-1.5-flash') # 2.5だと20回で制限にかかってしまうため
    
    prompt = f"""
    以下の論文のアブストラクトを読み、研究者向けに日本語で要約してください。
    
    【制約】
    ・3つの箇条書き（・を使用）で出力すること。
    ・専門用語はなるべくそのまま使用すること。
    ・結論を明確にすること。

    またこの論文の新規性や重要なポイントがあれば、要約の最後に追加してください。
    
    Abstract:
    {text}
    """
    
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        return f"要約生成エラー: {e}"