# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Documents and Settings\Administrator\.spyder2\.temp.py
"""

#!env python

import os
import pymssql
from math import * #导入数学运算，阶乘等
from decimal import Decimal
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties #font control
import numpy as np #np.roots()可以求高次方程的根
from scipy import stats
import cPickle as p
import time

#从cellCapacity表中获取小区的TCH信道数配置，经过calTradfficByTCH()计算，求出该小区
#信道配置情况下能支持的最大话务量erl，并写入临时表cellCapacityByTrafficTemp
#写入临时表是为了避免直接用PYTHON逐行update更行SQL表格，转而在SQL server中使用批量操作
target = 0.7         #usage to be approximate   

def data_drawer(cellTchs,cellP):
    conn = pymssql.connect(host="localhost",user="sql2008",password="123", database="ResourceAllocation", charset="utf8")
    cur = conn.cursor()

    #获取小区号和对应的P数目
    sql =  "select CELLID, busyVoiceTraffic from cellUsageByDayMax_Mini order by CELLID,time"
    cur.execute(sql)
    rows = cur.fetchall()
    allRecords = np.array(rows)#小区：小区中的tch数目

    for i in range(len(allRecords)):
        cellTchs[allRecords[i,0]] = allRecords[i,1]
    
    
    #获取小区号和对应的TCH
    sql =  "select CELLID, TCH from cellCapacity order by CELLID"
    cur.execute(sql)
    rows = cur.fetchall()
    allRecords = np.array(rows)#小区：小区中的tch数目

    for i in range(len(allRecords)):
        cellP[allRecords[i,0]] = allRecords[i,1]
    
    cur.close()
    conn.close() #close the connection to ms sql   

def get_init(P,Tchs,x,l):#produce the optimum point in Real Number Space
    #fileHandle = open ('test1.txt', 'w' )      
    
    c1=0.0
    c2=0.0
    for i in range(l):
        c1+=P[i]/target-Tchs[i]
        c2+=P[i]*P[i]
    c1*=8.0
    lamda=c1/c2
    print c1,c2,lamda
    for i in range(l):
        x[i]=-1*P[i]*P[i]*lamda/64.0+(P[i]/target-Tchs[i])/8.0
        #fileHandle.write(str(P[i])+' '+str(Tchs[i])+' '+str(x[i])+'\n')
    #fileHandle.close()

def f(x,l):
    sum = 0.0
    t = 0.0
    for i in range(l):
        t = x[i]-(P[i]/target-Tchs[i])/8
        sum+=64*t*t/(P[i]*P[i])
    return sum

#def descending_iteration(x,l):
#    y = x
#    for i in range(l):
        
def partial_f(x,P,Tchs):
    return 128.0*x/(P*P)-16/P/target+16*Tchs/(P*P)
    
def readDataFromFile(filename):
    fr = open(filename, 'r')
    data = p.load(fr)
    fr.close()
    return data


if __name__ == "__main__":
    start = time.clock()
    print "Time started: ", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

    #fileHandle = open ('test.txt', 'w' )  
    #writeCapacityToDB() #根据TCH配置数，计算个小区的容量，存入数据库临时表中
    cellTchs = {} # {小区：小区中的p数目}字典    
    cellP = {}    
    cellTchConf1 = readDataFromFile(r'cellTrafficByDayMax--LastMonth.data') 
    
    #print cellIDs
    
    cellTchConf2 = readDataFromFile(r'cellTchConfiguration.data') 
    cellIDs=cellTchConf2.keys()
    cellIDs.sort()
    print len(cellTchConf2)
    
    filename = "Tch1"
    fileHandle = open ("./data/"+filename+".txt", 'w')
    for id in cellIDs:
            fileHandle.write(str(int(id))+' '+str(cellTchConf2[id])+'\n')
#    for i in range(31):
#        filename = str(i+1)
#        fileHandle = open ("./data/P"+filename+".txt", 'w')
#        for id in cellIDs:
#            fileHandle.write(str(int(id))+' '+str(cellTchConf1[id][i])+'\n')
        
    #data_drawer(cellTchs,cellP)
#    Tchs = np.zeros(len(cellTchs),dtype=float)
#    P = np.zeros(len(cellTchs),dtype=float)
#    i = 0    
#    for key in cellTchs:
#        Tchs[i] = cellTchs[key]
#        P[i] = cellP[key]
#        i+=1
#    
#    l = len(P)
#    x = np.zeros(l,dtype=float)
#    get_init(P,Tchs,x,l)
#    y = np.zeros(l)
#    sum=0
#    for i in range(l):
#        #fileHandle.write(str(x[i])+'\n')
#        y[i] = round(x[i])
#        sum+=y[i]
#    #fileHandle.close()
#    print sum
    
    end = time.clock()
    print "Time ended: ",  time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())
    print "Time used: ", end-start, " seconds"
