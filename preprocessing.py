'''Prepares data for DL part of DeepXSOZ'''
import pandas as pd
import os, os.path
from tensorflow.keras.preprocessing.image import img_to_array
from numpy import save
from numpy import asarray
from numpy import load
import glob
from os import makedirs
from os import listdir
import shutil
from tensorflow.keras.preprocessing.image import load_img
import matplotlib.pyplot as plt
from os import path
import pathlib
import re
import tensorflow as tf
import time

class PreprocessData():
     """Extracts IC threshold images' labels of all subjects, and stores them in a
     folder "Subjects_NN_data" by relabelling them as Noise or Network"""
    
     def __init__(self, folder,subjectName):
        self.folder = folder
        self.subject = subjectName
        self.labels = self.create_label_list()
        self.relabel_IC_img()

    

     def create_label_list(self):
         """extract 0-1 labels and save them"""
        
         IC_labels = []
         folder_labels = self.folder
         

         #to sort files names numerically
         numbers = re.compile(r'(\d+)')
         def numericalSort(value):
             parts = numbers.split(value)
             parts[1::2] = map(int, parts[1::2])
             return parts
         
         #extract Noise (0) and Network (1,2,3) labels, convert labels > 1 as 1
         path = folder_labels + 'labels/' + self.subject
         print(path)

         for subject_label in sorted(glob.glob(path),key=numericalSort):
           

            print(subject_label)
            df = pd.read_csv(subject_label)

            #Label 0 - Noise, 1 - Normal RSN, 2 - Abnormal RSN, 3 - Atypical Signal Source
            df.loc[df['Label'] > 1, 'Label'] = 1
            labels = df.values[:, -1]
            print(labels)
            IC_labels.append(labels)
            flat_list = [item for sublist in IC_labels for item in sublist]
            
            
         #save all the labels 
         labels = asarray(flat_list)
         print(labels.shape)
         #save('NoiseNetworkLabels.npy', labels)
         return flat_list

     def relabel_IC_img(self):
        """ Relabel IC images of all the subjects as Noise (0)/Network (1) in a folder named 'Subjects_NN_data'"""

        #to sort files names numerically
        numbers = re.compile(r'(\d+)')
        def numericalSort(value):
            parts = numbers.split(value)
            parts[1::2] = map(int, parts[1::2])
            return parts

        #create a new folder 
        dataset = 'Subjects_NN_data/'
        makedirs(dataset, exist_ok=True)
        
        subjects = self.folder +  self.subject
        
        flag = False
      
        last_subject_ic_img_number = 0
        for subject in sorted(glob.glob(subjects),key=numericalSort):
              
             
            print(subject)
            subject_name = os.path.basename(os.path.normpath(subject))
            print("subjectname",subject_name)
            filename = subject
            num_ic_images_per_subject = len(glob.glob(filename + '/MO/report/f*.png'))
            print("num of images per sub",num_ic_images_per_subject)
            first_index = last_subject_ic_img_number
            last_index = num_ic_images_per_subject + first_index

            if(flag == True):
                first_index -= 1
                last_index = num_ic_images_per_subject + first_index
            subject_specific_label = self.labels[first_index:last_index]

            for ic_img in range(1,num_ic_images_per_subject+1):
         
                image_name = filename + '/MO/report/IC_' + str(ic_img) + '_thresh.png'
                shutil.copy(image_name, dataset)

                print(ic_img)
                if(subject_specific_label[ic_img-1] == 1):
                    destination = dataset + 'Network' + str(subject_name) + '_'+ str(ic_img)+ '_.png'
                    source = dataset + 'IC_' + str(ic_img) + '_thresh.png'
                    os.rename(source,destination)

                if(subject_specific_label[ic_img-1] == 0):
                    destination = dataset + 'Noise' + str(subject_name) + '_'+ str(ic_img)+ '_.png'
                    source = dataset + 'IC_' + str(ic_img) + '_thresh.png'
                    os.rename(source,destination)

            last_subject_ic_img_number = last_index +1
            flag = True
            

            
