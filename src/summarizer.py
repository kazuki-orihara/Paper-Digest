import google.generativeai as genai
from google.generativeai.types import HarmCategory, HarmBlockThreshold
import os

def summarize_abstract(abstract_text):
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        return "Error: Gemini API Key is missing."
        
    genai.configure(api_key=api_key)
    
    # Use a stable model version
    model = genai.GenerativeModel('gemini-2.5-flash')

    # Safety settings to prevent false positives on scientific content
    # (e.g., discussions on viruses, toxicity, or cellular death)
    safety_settings = {
        HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
    }

    # English Prompt for Global Audience
    # Changed 【 】 to ### for better Markdown compatibility
    prompt = f"""
    Please summarize the following research abstract for a professional biologist or researcher.
    Output the summary in English.
    
    ### Requirements
    1. **Background & Purpose**: Why was this study conducted?
    2. **Key Methodology**: What were the notable techniques or approaches?
    3. **Main Results**: What were the findings? (Include specific values or molecule names if available)
    4. **Conclusion & Impact**: What is the significance of this study to the field?
    
    ### Abstract
    {abstract_text}
    """
    
    try:
        response = model.generate_content(
            prompt,
            safety_settings=safety_settings
        )
        return response.text
    except Exception as e:
        # Return the error string so main.py can detect "429" or other issues
        return f"Summary Generation Error: {e}"