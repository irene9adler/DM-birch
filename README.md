# DM-birch

UCAS 刘莹-数据挖掘大作业 
birch算法c++实现
代码文件：1116_birch_debug  编译：sh gcc.sh   运行sh run.sh（可修改run.sh文件4个参数。）
run.sh参数设置：
-r + 文件名  -n + 数据条目数  -a + 数据属性数 -d + 半径阈值
示例：./birch -r 1.txt -n 1000 -a 2 -d 4
运行结果文件：out.txt 包括分类数目，各个类包含的items，运行时间。
