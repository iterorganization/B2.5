"""
Generates b2cdcn.F
"""

from __future__ import print_function

from xml.etree import ElementTree as ET
from xml.etree.ElementTree import tostring
import textwrap
from itertools import chain
import re
import os
import sys

def get_default(text):
    base = ' Defaults to '
    end = '.'
    if text:
        if '.true.' in text or '.false.' in 'nx' == text in 'ny' == text or \
            text == 'ns' in text:
            return base + text + end
        if text.isalpha() or text == ' ' or text.replace('.', '').isalpha() or\
           text.replace('_', '').isalpha():
            return base + "'" + text + "'" + end
        else:
            return base + text + end
    else:
        return ''

def replace_tags(text):
    """To be used before loading the XML text into ElementTree. BE CAREFUL ON
    WHAT TAGS YOU WISH TO REPLACE!!!

    Resolve <sup>*</sup>, <sub>*</sub>, ... tags before proceeding with
    parsing the XML:

    x<sup>a<sup> -> x^a
    x<sup>a,b,c<sup> -> x^{a,b,c}

    Arguments:
        text (str): Input text that contains tags we wish to replace.
    Returns:
        out (str): Output text with replaced tags.

    """
    # Define type of tags and their replacement. LaTeX notation.
    tags = {"sup": "^", "sub": "_"}
    out = text
    for tag in tags:
        pattern = f"<{tag}>(.+?)</{tag}>"
        matches = set(re.findall(pattern, text))
        for match in matches:
            replacement = f"{tags[tag]}"
            if len(match) > 1:
                # Add curly braces around. When using f-strings {{}} escapes
                # the curly braces.
                replacement += f"{{{match}}}"
            else:
                replacement += f"{match}"

            out = out.replace(f"<{tag}>{match}</{tag}>", replacement)
    return out

def dedent(description, prefix = '\n*    ', width=75):
    """ Removes first empty fort from description and any leading tabs
        from the next fort before the description and any following forts.
        First forts are wrapped to 70 characters.

    :param description(string): from the XML generated tooltips dictionary
    :return: formatted output for the tooltip
    """
    trim_start = 0  # Remove any leading newfort that affects dedent
    while trim_start < len(description) and description[trim_start] == '\n':
        trim_start += 1
    description = textwrap.dedent(description[trim_start:])
    lines = description.splitlines()
    output = ''
    for line in lines:
        line = line.lstrip()
        output += prefix[1:] + prefix.join(textwrap.wrap(line, width)) + '\n'
    return output[0:]  # remove last newline

def dedent_without_wraping(description):
    trim_start = 0  # Remove any leading newfort that affects dedent
    while trim_start < len(description) and description[trim_start] == '\n':
        trim_start += 1
    description = textwrap.dedent(description[trim_start:])
    return description

def add_to_fort(text):
    global fort
    split_text = text.splitlines()
    text = ''
    for i, line in enumerate(split_text):
        line = line.lstrip()
        if i==0:
            text += '\n*  ' + line[:2]
            text += '\n*     '.join(textwrap.wrap(line[2:], 74)) + '\n'
        else:
            text += '*     ' + '\n*     '.join(textwrap.wrap(line, 73)) + '\n'
    fort += text[:-1]

xml_name = "b2input.xml"

f = open(xml_name, 'r')
xml_text = f.read()
f.close()

"""PRE-PARSING

Resolve super-scripts and sub-scripts, which resides in <description> tags.
"""
xml_text = replace_tags(xml_text)

"""Controlling entities

Instead of writing greek letters as UTF-8 symbols for which pdflatex has
issues, generate a separate .F file that contains the latex commands for greek
letters.

Instead of hardcoding which are inside, use the:

set(re.findall("&(.+?);", text))

to get all unique commands.
"""

# The following flag deactivates resolving entities, namely Greek and math
# symbols from the definitions dtd file. Default entities such as &lt; &gt;
# &amp; are resolved by ElementTree
if 'DO_NOT_RESOLVE_ENTITIES' in os.environ:
    DO_NOT_RESOLVE_ENTITIES=1
else:
    DO_NOT_RESOLVE_ENTITIES=0

if DO_NOT_RESOLVE_ENTITIES:
    commands = set(re.findall('&(.+?);', xml_text))
    text_to_process = xml_text
    for command in commands:
        # The following entities are supported by default.
        if command in ["lt", "gt", "amp"]:
            continue
        text_to_process = text_to_process.replace(f'&{command};', f'\{command}')
else:
    dtd_name = "xhtml-symbol.ent"
    xml_entities = '<!ENTITY % symbols SYSTEM "xhtml-symbol.ent" > %symbols;'

    f = open(dtd_name, 'r')
    dtd_text = f.read()
    f.close()
    text_to_process = xml_text.replace(xml_entities, dtd_text)

tree = ET.ElementTree(ET.fromstring(text_to_process))
root = tree.getroot()

specification = root.find('module[@name="b2.parameters"]/routine[@name="b2cdcn"]')
purpose = specification.findtext('purpose')

DELIMITER = """*---------------------------------------------------------------\
--------\n"""

fort = DELIMITER
fort += """* Generated by b2cdcn.py from b2input.xml. DO NOT EDIT THIS FILE!
*.specification

      subroutine b2cdcn ()
      implicit none

"""

fort += DELIMITER
fort += '*.documentation\n*\n'
fort += '*  1. purpose\n*\n'
fort += dedent(purpose, prefix='\n*     ', width=66)
fort += '*\n' + DELIMITER
fort += '*.text\n*\n'


b2_parameters = root.find('module[@name="b2.parameters"]')
categories = b2_parameters.findall('category')


# Sections
for category in categories:


    # Namelist name
    namelist = category.findtext('introduction/namelist')
    fort += '* NAMELIST /' + namelist + '/\n'

    # Section description
    fort += '*    Found in ' + category.attrib['name']
    if category.findtext('introduction/note'):
        fort += category.findtext('introduction/note')
    fort += '.\n'
    fort += dedent(category.findtext('introduction/description'),
                   prefix='\n*    ')[:-1]
    # Adding all parameters and switches

    for element in category:
        if element.tag == 'switch':
            text = ''
            text += element.findtext('name') + ' - '
            text += element.findtext('type') + '. '
            text += dedent_without_wraping(element.find('description').text)
            text += get_default(element.findtext('default'))
            add_to_fort(text)
        elif element.tag == 'switchgroup':
            N = len(element.findall('switch'))
            counter = 0
            for switch in element.findall('switch'):
                counter+=1
                text = ''
                text += switch.findtext('name') + ' - '
                text += switch.findtext('type') + '. '
                text += dedent_without_wraping(switch.find('description').text)
                default = get_default(switch.findtext('default'))
                if default:
                    text += get_default(switch.findtext('default'))
                if counter == N and element.findtext('description'):
                    text += '\n'
                    text += dedent_without_wraping(element.find('description').text)

                add_to_fort(text)


        # Adding a standalone note.
        elif element.tag == 'note':
            fort += '*\n'
            text = dedent(element.text, prefix='\n* ')+'*\n'
            fort += text
            continue
        else:
            continue
    fort += '\n*\n'

fort += DELIMITER
fort += """*.end b2cdcn\n\n      end subroutine b2cdcn\n"""

if sys.version_info[0] >= 3:
    f = open('b2cdcn.F', 'wb')
    f.write(fort.encode())
else:
    f = open('b2cdcn.F', 'w')
    f.write(fort.encode('utf-8'))
f.close()
