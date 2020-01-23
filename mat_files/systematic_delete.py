"""
Delete files sytematically.
"""

import os
filename = input("filename: ")
#jump = int(input("jump: "))
#offset = int(input("offset: "))
#divisor = int(input("mod: "))
#extension = input("extension: ")

jump = 1000
offset = 1000
divisor = 4000
extension = "_1_1_HYmob.mat"

number = offset
while True:
    name = "%s%d%s"%(filename, number, extension)
    #"%s%d_0_1_HYmob.%s"%(filename, number, extension)
    if number % divisor == 0:
        os.remove(name)
        # print(name)
    number += jump
