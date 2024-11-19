# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:43:49 2023

@author: PetraH
"""

import cv2
import os
import config

if __name__ == "__main__":
    path_folder = '...' # path to folder with png files to convert to mp4
    
    image_folder = os.path.join(path_folder,'figs')
    video_name   = os.path.join(path_folder,'video.mp4')
    
    images = ['fig_t' + str(t)+'.png' for t in range(0,151)] # list figure names
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    
    video = cv2.VideoWriter(video_name, 0, 1, (width,height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()