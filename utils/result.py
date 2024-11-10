# https://chatgpt.com/c/672b9133-9d20-800e-894f-bd0889252edf

import os
from sympy import latex

class Result:
    def __init__(self, file_name="result.md", script_name="script.py", dir_path=""):
        self.config = {
            "file_name": file_name,
            "script_name": script_name,
            "dir_path": dir_path
        }

    def file_path(self):
        return os.path.join(self.config["dir_path"], self.config["file_name"])

    def create_result(self, mode="w"):
        # Create the file in 'x' mode (fails if file already exists)
        try:
            with open(self.file_path(), mode=mode) as file:
                file.write(f"*Executed from: {self.config['script_name']}*\n")
        except FileExistsError:
            print(f"File {self.file_path()} already exists.")
        return None

    def print_md(self, _string, end="\n"):
        """
        Save a result to the Markdown file.
        Appends the content to the file.
        """
        with open(self.file_path(), mode="a", encoding="utf-8") as file:
            file.write(_string + end)
        return None

    def print_latex(self, sympy_string):
        """
        Save a LaTeX formatted string to the Markdown file.
        """
        latex_string = f"${latex(sympy_string)}$\n"
        with open(self.file_path(), mode="a", encoding="utf-8") as file:
            file.write(latex_string)
        return None
