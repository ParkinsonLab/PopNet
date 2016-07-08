'''
Created on Jun 28, 2016

@author: javi
'''
import re

'''
Takes path
Returns (header, body) as separated lines
No further modifications
'''

def readTable(filepath):
    with open(filepath, 'r') as input:
        raw_text = input.read()
    
    line_split = re.split('\n', raw_text)
    header = line_split[0]
    body = line_split[1:-1]
    
    return header, body
        
'''
requires the table tuple (header, body) from readTable
returns the modified tuple in the same format
'''
def reverseOrder(table):        
    header, body = table
    
    grad_pattern = '(?<=colorlist=\"\")(.+?)(?=\")'
    val_pattern = '(?<=\")([0-9\.]+?[|].+?)(?=\")'
    results = []
    
    for line in body:
        grad = re.search(grad_pattern, line).group(1)
        vals = re.search(val_pattern, line).group(1)
        grad_list = re.split(',', grad)
        val_list = re.split('\|', vals)
        
        grad_list.reverse()
        val_list.reverse()
        
        new_grad = ','.join(grad_list)
        new_val = '|'.join(val_list)
        
        #sub the old strings with the new ones.
        mod_line = re.sub(grad_pattern, new_grad, line)
        mod_line2 = re.sub(val_pattern, new_val, mod_line)
        
        results.append(mod_line2)
    
    return (header, results)

def changeColor(input_path, color_dict):
    
    def colrep(matchobj):
        return color_dict[matchobj.group(0)]
    
    with open(input_path, 'r') as input:
        data = input.read()
    
    for color in color_dict:
        data = re.sub(color, colrep, data)
    
    return data
        
##helpers


if __name__ == '__main__':

    ## To Reverse Order
    directory = '/data/new/javi/reshape'
    file_name = 'yeast.csv'
    output_name = 'rev_yeast.csv'
    
    path = '/'.join([directory, file_name])
    output_path = '/'.join([directory, output_name])
    
    with open(output_path, 'w') as output:
        header, body = reverseOrder(readTable(path))
        output.write('\n'.join([header] + body))
        
#     ## To change color
#     directory = '/data/new/javi/reshape'
#     file_name = 'rev_toxo.csv'
#     output_name = 'rev_toxo_rc.csv'
#     
#     path = '/'.join([directory, file_name])
#     output_path = '/'.join([directory, output_name])
#     
#     color_dict = {'#FF9900' : '#0000FF',
#                   '#32FF00' : '#FF0000',
#                   '#3200FF' : '#00FF00',
#                   '#00FEFF' : '#FFa500'}
#     
#     with open(output_path, 'w') as output:
#         output.write(changeColor(path, color_dict))
    