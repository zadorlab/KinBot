"""CFOUR ZMAT formatting shared by generated and user-supplied inputs."""

from __future__ import annotations

import re


_KEYWORD_START = re.compile(r'(?im)^([ \t]*)\*CFOUR\(')
_MAX_KEYWORD_LINE = 72  # CFOUR 2.1 read only the first 80 columns in the live test.


def _split_keywords(body):
    """Split only top-level commas/newlines (e.g. CCSD(T) stays intact)."""
    words = []
    token = []
    depth = 0
    quote = None
    for char in body:
        if quote:
            token.append(char)
            if char == quote:
                quote = None
        elif char in ('"', "'"):
            quote = char
            token.append(char)
        elif char == '(':
            depth += 1
            token.append(char)
        elif char == ')':
            if depth == 0:
                raise ValueError('CFOUR keyword section has unmatched parentheses.')
            depth -= 1
            token.append(char)
        elif depth == 0 and char in ',\n':
            word = ''.join(token).strip()
            if word:
                words.append(word)
            token = []
        else:
            token.append(char)
    word = ''.join(token).strip()
    if word:
        words.append(word)
    if depth or quote or not words:
        raise ValueError('CFOUR keyword section is incomplete.')
    return words


def normalize_cfour_zmat(text):
    """Wrap *CFOUR keywords at safe line boundaries without changing values.

    CFOUR's manual says a newline separates keywords and forbids a trailing
    comma at the end of a continued line.  Each keyword gets its own line;
    a single unusually long keyword fails before submission.
    """
    output = []
    position = 0
    while match := _KEYWORD_START.search(text, position):
        output.append(text[position:match.start()])
        indent = match.group(1)
        begin = match.end()
        depth = 1
        quote = None
        end = None
        for index in range(begin, len(text)):
            char = text[index]
            if quote:
                if char == quote:
                    quote = None
            elif char in ('"', "'"):
                quote = char
            elif char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
                if depth == 0:
                    end = index
                    break
        if end is None:
            raise ValueError('CFOUR keyword section has no closing parenthesis.')
        keywords = _split_keywords(text[begin:end])
        lines = [indent + '*CFOUR(' + keywords[0]]
        lines.extend(indent + keyword for keyword in keywords[1:])
        lines[-1] += ')'
        if any(len(line) > _MAX_KEYWORD_LINE for line in lines):
            raise ValueError('CFOUR keyword exceeds the safe 72-column limit.')
        output.append('\n'.join(lines))
        position = end + 1
    output.append(text[position:])
    return ''.join(output)
