import os
import re
import argparse

def modify_link(line):
    replaced = False

    # Define icons and attributes
    icon_external_link = ':material-open-in-new:'
    icon_download = ':material-download:'
    open_new_tab = '{:target="_blank"}'

    # Define patterns for general and download links
    general_pattern = r'\[([^\]]+)\]\((http[s]?://[^\s\)]+)\)'
    download_pattern = r'\[([^\]]+)\]\((https://polybox\.ethz\.ch/index\.php/s/[^\s\)]+)\)'

    # Define replacements for general and download links
    general_replacement = r'[\1 ' + icon_external_link + r'](\2)'
    download_replacement = r'[\1 ' + icon_download + r'](\2)'

    # Check if the general pattern appears in the line
    if re.search(general_pattern, line):
        # Check for link icon
        if icon_external_link not in line and icon_download not in line:
            new_line = re.sub(download_pattern, download_replacement, line)
            if new_line != line:
                line = new_line
                replaced = True
            else:
                new_line = re.sub(general_pattern, general_replacement, line)
                if new_line != line:
                    line = new_line
                    replaced = True

        # Check for new tab attribute
        if open_new_tab not in line:
            new_line = re.sub(r'(\[.*?\]\(.*?\))', r'\1' + open_new_tab, line)
            if new_line != line:
                line = new_line
                replaced = True

    return line, replaced


def process_markdown_file(file_path):
    with open(file_path, 'r+', encoding='utf-8') as file:
        lines = file.readlines()
        modified = False
        for i, line in enumerate(lines):
            new_line, changed = modify_link(line)
            if changed:
                modified = True
                lines[i] = new_line
                print(f"Modifying line {i+1} in file: {file_path}")
        if modified:
            file.seek(0)
            file.writelines(lines)
            file.truncate()
            print(f"File modified: {file_path}")


def main(start_path):
    for root, dirs, files in os.walk(start_path):
        for file in files:
            if file.endswith('.md'):
                process_markdown_file(os.path.join(root, file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check markdown links in files.')
    parser.add_argument('-p', '--path', default=os.getcwd(), help='Base path to search for markdown files. Defaults to current working directory.')
    args = parser.parse_args()

    main(args.path) 
