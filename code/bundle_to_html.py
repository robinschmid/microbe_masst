from bs4 import BeautifulSoup
from pathlib import Path
import base64
import requests
import sys
import argparse
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def replace_by_local_file(path):
    """
    Replaces web ressources by local files if possible
    :param path: URL
    :return: a local file or the original input string
    """
    if path == "https://d3js.org/d3.v3.min.js":
        return "../code/lib/d3.v3.min.js"
    if path == "https://d3js.org/d3.v6.min.js":
        return "../code/lib/d3.v6.min.js"
    if path == "https://code.jquery.com/jquery-3.6.0.min.js":
        return "../code/lib/jquery-3.6.0.min.js"
    if str(path).endswith("jquery.dataTables.min.js"):
        return "../code/lib/jquery.dataTables.min.js"
    if str(path).endswith("jquery.dataTables.min.css"):
        return "../code/lib/jquery.dataTables.min.css"
    return path


def replace_data_in_file(data_json_file, text, placeholder_str):
    """

    :param data_json_file: data might be passed as file or other data structure like json
    :param text: the text to search and replace in
    :param placeholder_str: the placeholder string to replace
    :return: the final text
    """
    try:
        tree_data = Path(data_json_file).read_text()
    except:
        tree_data = data_json_file
        pass
    text = text.replace(placeholder_str, tree_data, 1)

    return text


def build_dist_html(input_html, output_html, replace_dict=None, compress=False):
    """
    Creates a single distributable HTML file.
    Reads the input_html and internalizes all CSS, JS, and data files into the output html. For web ressources: First
    try to load a local file, else try to download file.
    :param replace_dict: dict of key (placeholder string in HTML) and the value to insert either from a file or as a
    string
    :param input_html: the input html file that defines all dependencies
    :param output_html: the bundled HTML file
    :param data_json_file: data as a file or as json (or dictionary / array)
    :param placeholder_str: the placeholder is replaced with the data
    :return: None
    """
    if replace_dict is None:
        replace_dict = {}
    original_html_text = Path(input_html).read_text(encoding="utf-8")
    soup = BeautifulSoup(original_html_text, "lxml")

    # Find link tags. example: <link rel="stylesheet" href="css/somestyle.css">
    for tag in soup.find_all('link', href=True):
        if tag.has_attr('href'):
            path = tag['href']
            path = replace_by_local_file(path)
            file_text = Path(path).read_text(encoding="utf-8")

            # remove the tag from soup
            tag.extract()

            # insert style element
            new_style = soup.new_tag('style')
            new_style.string = file_text
            soup.html.head.append(new_style)

    # Find script tags. example: <script src="js/somescript.js"></script>
    for tag in soup.find_all('script', src=True):
        if tag.has_attr('src'):
            path = tag['src']
            path = replace_by_local_file(path)
            if path.startswith("http"):
                response = requests.get(path)
                response.raise_for_status()
                file_text = response.text
            else:
                file_text = Path(path).read_text()

            # remove the tag from soup
            tag.extract()

            # insert script element
            new_script = soup.new_tag('script')
            new_script.string = file_text
            soup.body.append(new_script)

    # Find image tags.
    for tag in soup.find_all('img', src=True):
        if tag.has_attr('src'):
            file_content = Path(tag['src']).read_bytes()

            # replace filename with base64 of the content of the file
            base64_file_content = base64.b64encode(file_content)
            tag['src'] = "data:image/png;base64, {}".format(base64_file_content.decode('ascii'))

    out_text = str(soup)

    # try to replace data with PLACEHOLDER_JSON_DATA
    for placeholder, replace_data in replace_dict.items():
        if replace_data is not None:
            out_text = replace_data_in_file(replace_data, out_text, placeholder)
        else:
            out_text = replace_data_in_file("", out_text, placeholder)

    if compress:
        try:
            import minify_html
            out_text = minify_html.minify(out_text, minify_js=True, minify_css=True)

        except Exception as e:
            logger.warning("Error during output compression.")
            logger.exception(e)

    # Save onefile

    with open(output_html, "w", encoding="utf-8") as outfile:
        outfile.write(out_text)

    return True


if __name__ == '__main__':
    # parsing the arguments (all optional)
    parser = argparse.ArgumentParser(description='Parse an HTML file and all dependencies to create a single '
                                                 'distributable HTML file')
    parser.add_argument('--in_html', type=str, help='The input html file',
                        default="../code/collapsible_tree_v3.html")
    parser.add_argument('--replace_dict', type=str, help='Replace placeholder strings with values or file contents. '
                                                         'form of {"PLACEHOLDER": "FILE OR STRING"}',
                        default="../output/merged_ontology_data.json")
    parser.add_argument('--out_html', type=str, help='output html file', default="../output/oneindex.html")
    parser.add_argument('--compress', type=bool, help='Compress output file (needs minify_html)',
                        default=True)
    args = parser.parse_args()

    # is a url - try to download file
    # something like https://raw.githubusercontent.com/robinschmid/GFOPontology/master/data/GFOP.owl
    # important use raw file on github!
    try:
        build_dist_html(args.in_html, args.out_html, args.replace_dict, args.compress)
    except Exception as e:
        # exit with error
        logger.exception(e)
        sys.exit(1)

    # exit with OK
    sys.exit(0)
