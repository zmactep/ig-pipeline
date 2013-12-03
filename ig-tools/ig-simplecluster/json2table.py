from argparse import ArgumentParser
from functools import partial
import json
from xml.etree.ElementTree import TreeBuilder, ElementTree
import sys

class JsonKeyError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('json_filename', help='input json filename')
    parser.add_argument('headers_filename', help='table column headers json filename')
    parser.add_argument('out_filename', help='output html filename')
    return parser.parse_args()


class Node:
    def __init__(self, text, *children):
        self.children = children
        self.text = text
        self.is_leaf = not children  # i.e. empty
        self.width = self.measure()

    def measure(self):
        children_width = 0
        for c in self.children:
            c.measure()
            children_width += c.width
        return max(1, children_width)


def render_data(builder, trees):
    def render(node, row_opened):
        if not row_opened:
            builder.start('tr')
        builder.start('td', {'rowspan': str(node.width)})
        builder.data(node.text)
        builder.end('td')
        for c in node.children[:1]:
            render(c, True)
        for c in node.children[1:]:
            render(c, False)
        if node.is_leaf:
            builder.end('tr')
    for root in trees:
        render(root, False)


def render_headers(builder, trees):
    def render(nodes):
        if nodes:
            builder.start('tr')
            for node in nodes:
                builder.start('th', {'colspan': str(node.width)})
                builder.data(node.text)
                builder.end('th')
            builder.end('tr')
            render([child for node in nodes for child in node.children])
    render(trees)


def build_forest(data):
    forest = []
    if type(data) == dict:
        for k, v in data.items():
            forest.append(Node(k, *build_forest(v)))
    elif type(data) == list:
        for v in data:
            forest.append(*build_forest(v))
    else:
        assert type(data) == str
        forest.append(Node(data))

    return forest


def load_json_value(json_filename, key):
    with open(json_filename) as file:
        data = json.load(file)
    if not key in data:
        raise JsonKeyError('No \'{0}\' root element key in {1}.'.format(key, json_filename))
    return data[key]


def load_headers(json_filename):
    return build_forest(load_json_value(json_filename, 'labels'))


def load_data(json_filename):
    return build_forest(load_json_value(json_filename, 'groups'))


def get_xhtml_tree(builder, header_generator, row_generator, title):
    builder.start('html')
    builder.start('head')
    builder.start('title')
    builder.data(title),
    builder.end('title')
    builder.end('head')
    builder.start('body')
    builder.start('table', {'border': '1'})
    header_generator(builder)
    row_generator(builder)
    builder.end('table')
    builder.end('body')
    builder.end('html')
    return ElementTree(builder.close())


def main():
    args = parse_args()

    try:
        data_forest = load_data(args.json_filename)
        headers_forest = load_headers(args.headers_filename)
    except JsonKeyError as e:
        print(e.message)
        sys.exit(-1)

    if not data_forest:
        print('No matching root node in input data.')
        sys.exit(-1)

    builder = TreeBuilder()

    with open(args.out_filename, 'w') as f:
        get_xhtml_tree(builder,
                       partial(render_headers, trees=headers_forest),
                       partial(render_data, trees=data_forest),
                       '{0} as table'.format(args.json_filename))\
            .write(f, encoding="unicode")


if __name__ == '__main__':
    main()