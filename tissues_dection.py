from tkinter import W
import openslide as ops 
import cv2
import io 
import sys
import numpy as np
from PIL import Image
import multiprocessing
#from lab_nor import *
import pdb
import functools
import random
FLAG = 0


# def histeq(imarr):
#   hist, bins = np.histogram(imarr, 100)
#   cdf = np.cumsum(hist)
#   cdf = 100 * (cdf/cdf[-1])
#   res = np.interp(imarr.flatten(), bins[:-1], cdf)
#   res = res.reshape(imarr.shape)
#   return res, hist

def identify(axis_tup, slide_path, output_base, size, cut, thr):
    #通过RGB标准差值判断是否为组织
    
    w = axis_tup[0] * (size - cut)
    h = axis_tup[1] * (size - cut)
    #重叠cut
    slide = ops.open_slide(slide_path)
    region = slide.read_region((w,h), 0 , (size, size)).convert('RGB')
    region_array = np.array(region)
    
    # for i in range(size):
    #         for j in range(size):
    #             if((region_array[i,j,:] == [0,0,0]).all()):
    #                 region_array[i,j,:] = [255, 255, 255]

    var_R = np.std(region_array[:,:,0])
    var_G = np.std(region_array[:,:,1])
    var_B = np.std(region_array[:,:,2])
    # print(str((var_R + var_G + var_B)/3))

    if((var_R + var_G + var_B)/3 > thr):
    
        # w = region.size[0]
        # h = region.size[1]
        # region_new = np.zeros((w,h,3))
        # lab = np.zeros((w,h,3))
        # for i in range(w):
        #     for j in range(h):
        #         Lab = RGB2Lab(patch_array[i,j,:])
        #         lab[i, j] = (Lab[0], Lab[1], Lab[2])
        #均衡化
        #print("max_L:{:.2%},min_L:{:.2%}".format(max(lab[:,:,0]),min(lab[:,:,0])))
        #pdb.set_trace()
        #l_res, l_hist = histeq(lab[:,:,0])
        
        #print("max_eL:{:.2%},min_eL:{:.2%}".format(max(equ_L),min(equ_L)))
        #lab[:,:,0] = l_res
 
        # for i in range(w):
        #     for j in range(h):
        #         rgb = Lab2RGB(lab[i,j])
        #         region_new[i, j] = (rgb[0], rgb[1], rgb[2])
 
	    #cv2.imwrite(r'E:\code\collor_recorrect\test.jpg', img_new)
        # img = Image.fromarray(np.uint8(region_new)).convert('RGB')
        global FLAG
        FLAG += 1
        if FLAG<3000:
            print(FLAG)
            im = Image.fromarray(region_array)
            im.save(output_base+"/tile_x_"+str(w+20)+"_y_"+str(h+20)+"_szx_"+str(size-cut)+"_szy_"+str(size-cut)+".jpg")
        # print("save"+str((var_R + var_G + var_B)/3))

# def save_image(save_tup, slide_path, output_base, size, cut):

#     w = save_tup[0]
#     h = save_tup[1]
#     slide = ops.open_slide(slide_path)
    
#     region = slide.read_region((w,h), 0 , (size, size)).convert('RGB')
#     region_array = np.array(region)

#     im = Image.fromarray(region_array)
#     im.save(output_base+"/tile_x_"+str(w+20)+"_y_"+str(h+20)+"_szx_"+str(size-cut)+"_szy_"+str(size-cut)+".jpg")


def main():

    if len(sys.argv)<2:
        print("参数个数不足")
        sys.exit(-1)

    slide_path = sys.argv[1]
    output_base = sys.argv[2]
    # mask_path = sys.argv[3]
    
    # slide_path = "/data/liangjianwen/Mutil_GTEx/breast/gtex_breast_wsi/GTEX-1117F-2825.svs"
    # output_base = "/data/liangjianwen/Mutil_GTEx/breast/OUT_TILE/GTEX-1117F-2825-seg-parallel"
    # mask_path = "/data/liangjianwen/Mutil_GTEx/breast/OUT_MASK/GTEX-1117F-2825_mask.png"
    

    
    slide = ops.open_slide(slide_path)
    ww = slide.dimensions[0]
    hh = slide.dimensions[1]
    #slide = slide.read_region((0,0), 0 , (ww, hh)).convert('RGB')
    s = 500
    c = 40
    # slide_array = np.array(slide)
    # mask_5x  = Image.open(mask_path).convert('L')

    # mask_20x = mask_5x.resize((ww, hh), resample=(0))
    # mask_array = np.array(mask_20x, dtype=np.int)/255
    # slide_array[:,:,0] = slide_array[:,:,0] * mask_array
    # slide_array[:,:,1] = slide_array[:,:,1] * mask_array
    # slide_array[:,:,2] = slide_array[:,:,2] * mask_array
    # print(' mul done!')
    #print((set(mask_array.flatten().tolist())))
    # for tup in [(i,j) for i in range(int((ww-s)/(s-c))+1) for j in range(int((hh-s)/(s-c))+1)]:
    #     identify(axis_tup=tup, slide_array=slide_array, output_base=output_base, size=s, cut=c, thr=42)
    iter = [(i,j) for i in range(int((ww-s)/(s-c))+1) for j in range(int((hh-s)/(s-c))+1)]
    print(len(iter))
    index_list = []
    iter_3000 = []
   
    index_list = random.sample(range(len(iter)), len(iter))
    iter_3000 = [iter[i] for i in index_list]
    
    iter = iter_3000
    # save_list = []

    # ncore = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(processes=(60))
    # print("可用核数："+str(ncore))



    # par_identify = functools.partial(identify, slide_path=slide_path, output_base=output_base, size=s, cut=c, thr=22)
    # pool.map(par_identify, iter)
    # pool.close()
    # pool.join()

    for it in iter:
        identify(axis_tup=it,slide_path=slide_path, output_base=output_base, size=s, cut=c, thr=22)

    
    
   
    

if __name__ == '__main__':
    main()