orig文件夹：原始数据集图像
real-bin文件夹：ground truth
目前的result文件夹：我们的最后结果，就是用来计算的。

①edge_extraction.m:检测orig文件夹原始图像的边缘，生成的图片保存在orig-operation文件夹
②line_fitting.m:检测并去除orig-operation文件夹中图片的直线，生成的图片保存在remove-line文件夹
③ellipse_fitting.m:检测remove-line文件夹中的椭圆，并将检测出的椭圆保存在result425文件夹，
之后根据其与real-bin文件夹中对应的ground truth计算检测效果。③的调用函数包括④⑤⑥：
④fit_ellipse.m，⑤fit_ellipse_state.m，⑥fit_ellipse_state2.m:椭圆生成算法
多次运行③，将result425中得到的较好结果转移到'目前的result'文件夹
⑦shizuo.m:根据'目前的result'文件夹和orig文件夹的原始图像，将检测出来的椭圆边线加到原始图像上,
并膨胀。

曾昕阳，19.5.5
修改于19.6.18

