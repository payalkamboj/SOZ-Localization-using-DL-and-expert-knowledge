import numpy as np
from os import listdir
from os import makedirs
from os import listdir
from shutil import copyfile
from random import seed
from random import random
import sys
from matplotlib import pyplot
from keras.models import Sequential
from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import Flatten
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.optimizers import Adam
import tensorflow as tf
from keras import callbacks


class DLNoiseCLassification():
	
	def __init__(self, folder):
		 self.datasetHome = folder

	def CreateTrainTestDataset(self):
		dataset_home = self.datasetHome
		subdirs = ['train/', 'test/']
		for subdir in subdirs:
		# create label subdirectories
		#written Nqetwork instead of network to make it label 1 and Noise as label 0
			labeldirs = ['Nqetwork/', 'Noise/']
			for labldir in labeldirs:
				newdir = dataset_home + subdir + labldir
				makedirs(newdir, exist_ok=True)
		seed(1)
		# define ratio of pictures to use for validation
		val_ratio = 0.25
		# copy training dataset images into subdirectories
		src_directory = 'Subjects_NN_data/'
		for file in listdir(src_directory):
			src = src_directory + '/' + file
			dst_dir = 'train/'
			if random() < val_ratio:
				dst_dir = 'test/'
			if file.startswith('Noise'):
				dst = dataset_home + dst_dir + 'Noise/'  + file
				copyfile(src, dst)
			elif file.startswith('Network'):
				dst = dataset_home + dst_dir + 'Nqetwork/'  + file
				copyfile(src, dst)

	
	def TrainingTestingCNN(self):
		#define cnn model
		def define_model():
			model = Sequential()

			model.add(Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_uniform', padding='same', input_shape=(270, 400, 3)))
			model.add(MaxPooling2D((2, 2)))
			model.add(Conv2D(64, (3, 3), activation='relu', kernel_initializer='he_uniform', padding='same'))
			model.add(MaxPooling2D((2, 2)))
			model.add(Conv2D(256, (3, 3), activation='relu', kernel_initializer='he_uniform', padding='same'))
			model.add(MaxPooling2D((2, 2)))
			model.add(Flatten())
			model.add(Dense(704, activation='relu', kernel_initializer='he_uniform'))
			model.add(Dropout(0.33))
			model.add(Dense(1, activation='sigmoid'))

			# compile model
			opt = Adam(learning_rate=0.0001)
			model.compile(optimizer=opt, loss='binary_crossentropy', metrics=['accuracy', tf.keras.metrics.FalseNegatives(),tf.keras.metrics.FalsePositives(),tf.keras.metrics.TrueNegatives(),tf.keras.metrics.TruePositives()])
			return model


		# run the test harness for evaluating a model
		def run_test_harness():
			# define model
			model = define_model()
	
			# create data generator
			datagen = ImageDataGenerator(rescale=1.0/255.0,validation_split=0.1)
	
	
			# prepare iterators
			train_it = datagen.flow_from_directory('CNNDataTrainTest/train/',class_mode='binary', target_size=(270, 400))
			val_it = datagen.flow_from_directory('CNNDataTrainTest/train/',class_mode='binary', target_size=(270, 400),subset = 'validation')
			#shuffle=False to preserve order of files names
			test_it = datagen.flow_from_directory('CNNDataTrainTest/test/',class_mode='binary', target_size=(270, 400),shuffle=False)

			earlystopping = callbacks.EarlyStopping(monitor ="val_loss", 
												mode ="min", patience = 3, 
												restore_best_weights = True) 
	
			# fit model
	
			history = model.fit(train_it, epochs=20,validation_data = val_it, verbose=0,callbacks =[earlystopping])
			model.save('CNN_modelTrained')
	
			# evaluate model
			_, acc,FN, FP, TN, TP = model.evaluate(test_it, verbose=1)
			print("Accuracy: "'> %.3f' % (acc * 100.0))
			print("False_Negative\n", FN)
			print("False_Positive\n", FP)
			print("True_Negative\n", TN)
			print("True_Positive\n", TP)

			filenames = test_it.filenames
			nb_samples = len(filenames)
			predict = model.predict_generator(test_it,steps = nb_samples)
			y_pred = np.rint(predict)
			y_true = test_it.classes
			columns = []
			columns += [filenames,y_pred,y_true]
			np.savetxt("CNNNoiselabels.txt", np.column_stack(columns),fmt='%s', delimiter="  ")

		# entry point, run the test harness
		run_test_harness()
