from numpy import asarray
from numpy import save
from numpy import load
from tensorflow.keras.models import Sequential
from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from keras.layers import Dense
from tensorflow.keras.layers import Flatten
from keras.layers import Dropout
from keras.preprocessing.image import ImageDataGenerator
import tensorflow as tf
from tensorflow import keras
from os import listdir
import keras_tuner as kt


def build_model(hp):
    
   

    model = Sequential()
    
    model.add(Conv2D(filters=hp.Choice('Input_layer', values=[16,32,64,128,256, 512,1024,1500,2048], default=128,), kernel_size = (3, 3),input_shape = (270, 400, 3), activation = 'relu', name = "convolution_layer_1"))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Conv2D(filters=hp.Choice('num_filters1', values=[16,32,64,128,256, 512,1500,2048], default=128,), kernel_size = (3, 3), activation = 'relu', name = "convolution_layer_2"))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Conv2D(filters=hp.Choice('num_filters2', values=[16,32,64,128,256, 512,1500,2048], default=128,), kernel_size = (3, 3), activation = 'relu', name = "convolution_layer_3"))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Conv2D(filters=hp.Choice('num_filters3', values=[16,32,64,128,256, 512,1500,2048], default=128,), kernel_size = (3, 3), activation = 'relu', name = "convolution_layer_4"))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Conv2D(filters=hp.Choice('num_filters3', values=[16,32,64,128,256, 512,1500,2048], default=128,), kernel_size = (3, 3), activation = 'relu', name = "convolution_layer_4"))
    model.add(MaxPooling2D(pool_size = (2, 2)))
    model.add(Flatten())
    model.add(Dense(hp.Int("dense_units", min_value=192,max_value=4096, step=256),activation="relu"))
    model.add(Dropout(0.33))
    model.add(Dense(1, activation='sigmoid'))
    

    model.compile(optimizer=keras.optimizers.Adam(hp.Choice("learning_rate", values=[1e-2, 1e-3, 1e-4])),
        loss="binary_crossentropy",
        metrics=["accuracy"])
    return model

#Hyperband tuner
tuner = kt.Hyperband(
    build_model,
    objective='val_accuracy',
    max_epochs=20,project_name='CNN_tuner_highParameters')

datagen = ImageDataGenerator(rescale=1.0/255.0,validation_split=0.12)
	
# prepare iterators
train_it = datagen.flow_from_directory('CNNDataTrainTest/train/',class_mode='binary', target_size=(270, 400))
val_it = datagen.flow_from_directory('CNNDataTrainTest/train/',class_mode='binary', target_size=(270, 400),subset = 'validation')

earlystopping = tf.keras.callbacks.EarlyStopping(monitor ="val_loss", 
                                        mode ="min", patience = 5, 
                                        restore_best_weights = True) 
tuner.search(train_it, epochs=20, verbose=1,validation_data = val_it)
best_model = tuner.get_best_models()[0]
tuner.results_summary()