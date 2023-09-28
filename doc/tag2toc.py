import xml.etree.ElementTree as ET
from pathlib import Path
import sys
from collections import defaultdict
from textwrap import dedent


tagfile = "build/doc/sphinx/cpp/html/tagfile.xml"
tree = ET.parse(tagfile)
root = tree.getroot()

filenames = {}
group_classes = defaultdict(list)
group_subgroups = defaultdict(list)
group_titles = {}
child_groups = set()

for xgroup in root.findall("compound[@kind='group']"):
    name = xgroup.find("name").text
    filenames[name] = Path(xgroup.find("filename").text).stem
    group_titles[name] = xgroup.find("title").text
    for grp in xgroup.findall("subgroup"):
        group_subgroups[name].append(grp.text)
        child_groups.add(grp.text)
    for cls in xgroup.findall("class"):
        group_classes[name].append(cls.text)

for xgroup in root.findall("compound[@kind='class']"):
    name = xgroup.find("name").text
    filenames[name] = Path(xgroup.find("filename").text).stem

with open("build/doc/sphinx/cpp/html/modules.rst", "a") as out:
    out.write(dedent("""
              .. toctree::
                 :maxdepth: 2
                 :hidden:

              """))

    for group, title in group_titles.items():
        if group not in child_groups:
            out.write(f"   {title} <{filenames[group]}>\n")

for group, title in group_titles.items():
    with open(f"build/doc/sphinx/cpp/html/{filenames[group]}.rst", "a") as out:
        if any(subgroup in filenames for subgroup in group_subgroups[group]):
            out.write(dedent("""
                      .. toctree::
                         :caption: Groups
                         :maxdepth: 2
                         :hidden:

                    """))

            for subgroup in group_subgroups[group]:
                if subgroup in filenames:
                    out.write(f"   {group_titles[subgroup]} <{filenames[subgroup]}>\n")

        if any(cls in filenames for cls in group_classes[group]):
            out.write(dedent("""
                      .. toctree::
                         :caption: Classes
                         :maxdepth: 2
                         :hidden:

                    """))

            for cls in group_classes[group]:
                if cls in filenames:
                    out.write(f"   {cls[9:]} <{filenames[cls]}>\n")
