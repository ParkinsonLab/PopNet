'''
Created on Apr 29, 2014

@author: javi
'''

import random

def toTail():
    num = 0
    count = 0
    while num != 1:
        num = random.randint(0,400)
        count += 1
    return count

def pledge(num):
    count = 0. 
    for x in range(num):
        if random.randint(0,400) == 1:
            count += 1
    return count/num

if __name__ == '__main__':
    stuff = []
    reps = 100000
    for x in range(1,reps):
        stuff.append(toTail())
    print(stuff[0], stuff[len(stuff)-1])
    print("median:" + str(sorted(stuff)[len(stuff)/2]))
    
    average = 0. 
    reps = 100000
    for x in range(1, reps):
        average = (average * (x-1) + toTail()) / x
    print(average)