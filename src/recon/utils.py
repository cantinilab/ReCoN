import re

def split_layer_name(s, separator='_'):
    # Use re.escape to handle special characters in the separator
    pattern = rf"(.*){re.escape(separator)}(\w+)$"
    match = re.match(pattern, s)
    if match:
        return match.group(1), match.group(2)
    else:
        return s, ''  # In case there's no separator, return the original string and empty

