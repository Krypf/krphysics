# import os

def print_md(_string):
    """
    Save a result to 'result.md' in the root directory.
    """
    # Append to result.md (create the file if it doesn't exist)
    with open('result.md', 'w') as file:
        file.write(_string)
    return None

def print_latex(latex_string, script_name='script.py'):
    # Prepare the string to append
    latex_string = f"${latex_string}$\n"
    
    # Append to result.md (create the file if it doesn't exist)
    with open('result.md', 'a') as file:
        file.write(f"\n*Executed from: {script_name}*\n")
        file.write(latex_string)
    return None