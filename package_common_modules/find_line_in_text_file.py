'''
'''
# *** Find_line finds a line in .txt file which contains 'needed_text'  ***
#     needed_text - text to find
#     start_char - number of character where the needed_text starts in line
#     line - returns a text line from the file with needed_text
def find_line_in_text_file (file_handle, needed_text, start_char, line):
    tempChar = ''
    while tempChar != needed_text:
        line = file_handle.readline()
        tempChar = line[start_char : len(needed_text) + start_char]
    return line
