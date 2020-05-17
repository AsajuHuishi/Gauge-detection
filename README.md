## Gauge-detection
It is the raw code for Automatic Gauge Detection via Geometric Fitting for Safety Inspection. ([Li B, Yang J, Zeng X, et al. Automatic Gauge Detection via Geometric Fitting for Safety Inspection[J]. IEEE Access, 2019, 7: 87042-87048.](https://ieeexplore.ieee.org/document/8746263))

## Requirements
MATLAB

## Introduction
Folder:
orig：original image 
real-bin：ground truth
目前的result文件夹：our result for calculate

m File:
①edge_extraction.m:检测orig文件夹原始图像的边缘，生成的图片保存在orig-operation文件夹
②line_fitting.m:检测并去除orig-operation文件夹中图片的直线，生成的图片保存在remove-line文件夹
③ellipse_fitting.m:检测remove-line文件夹中的椭圆，并将检测出的椭圆保存在result425文件夹，
之后根据其与real-bin文件夹中对应的ground truth计算检测效果。③的调用函数包括④⑤⑥：
④fit_ellipse.m，⑤fit_ellipse_state.m，⑥fit_ellipse_state2.m:椭圆生成算法
多次运行③，将result425中得到的较好结果转移到'目前的result'文件夹
⑦shizuo.m:根据'目前的result'文件夹和orig文件夹的原始图像，将检测出来的椭圆边线加到原始图像上,并膨胀。

## Results
1.Comparison of detected shapes of gauges with different
methods. From left to right: the input gauge images, CHT method and our method:
<div align='center'>
<img src='https://img-blog.csdnimg.cn/20200517192922984.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzM2OTM3Njg0,size_16,color_FFFFFF,t_70' label='Comparison of detected shapes of gauges with different
methods. From left to right: the input gauge images, CHT method and our
method.'>
</div>

2.Steps of our gauge detection method: a) the input gauge image, b) result of edge detection, c) result of the line fitting, d) result of
line removal, e) result of ellipse fitting and d) the superposition of the result of ellipse fitting and the original input gauge image

<div align='center'>
<img src='https://img-blog.csdnimg.cn/20200517194108735.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzM2OTM3Njg0,size_16,color_FFFFFF,t_70'>
</div>

3.Visual comparison results. a) Ours, b) CHT, c) CCOEFF, d) CCOEFF_NORM, e) CCORR, f) CCORR_NORM, g) SQDIFF and
h) SQDIFF_NORM
<div align='center'>
<img src='https://img-blog.csdnimg.cn/20200517194531340.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3FxXzM2OTM3Njg0,size_16,color_FFFFFF,t_70'>
</div>