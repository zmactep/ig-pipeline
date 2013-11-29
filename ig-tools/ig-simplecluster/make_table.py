__author__ = 'pavel'

import json


def json2table(dt):
    d = json.loads(dt)["groups"]

    text = """
    <html>
    <body>
    <table border="1">
        <th><td colspan="2">CDR<td><td colspan="2">Full</td></th>
        <tr><td>Cluster<td><td>Similarity</td><td>Similarity</td><td>Clone</td></tr>
        {}
    </table>
    </body>
    </html>
    """

    to_append = []
    for cluster in d:
        scc = []
        for sim_cdr in d[cluster]:
            sfc = []
            for sim_full in d[cluster][sim_cdr]:
                lc = []
                for clone in d[cluster][sim_cdr][sim_full]:
                    lc.append("<td style=\"font-size: 8pt;\">%s</td><tr>" % clone)
                t = "".join(lc)
                sfc.append("<td colspan=\"{}\" style=\"font-size: 8pt;\">{}</td>".format(len(lc), sim_full) + t)
            t = "".join(sfc)
            scc.append("<td colspan=\"{}\" style=\"font-size: 8pt;\">{}</td>".format(len(sfc), sim_cdr) + t)
        t = "".join(scc)
        to_append.append("<tr><td colspan=\"{}\" style=\"font-size: 8pt;\">{}</td>".format(len(scc), cluster) + t)

    return text.format("".join(to_append))


def main():
    import sys
    js = open(sys.argv[1], "rt").readline()
    with open(sys.argv[2], "wt") as fd:
        fd.write(json2table(js))

if __name__ == "__main__":
    main()