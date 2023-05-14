# Finding overlapping SV calls
We have called and filtered SVs using three different tools! Now we take a sneak peak at how well the tools agree with one another. For this exercise, we're going to use [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR). This is a tool that considers SV type, size, location, strand and overall agreement among multiple SV calling tools. 

To get started, let's add SURVIVOR to our `$PATH`. *BE CAREFUL HERE*. We're adding something to our bash path to make it easy to call SURVIVOR, but you can overwrite your ability to call on other programmes. So be mindful when altering your `$PATH`. 
```
export PATH="${PATH}:/home/SURVIVOR/Debug/"
```
Now we should be able to call SURVIVOR with a simple:
```
SURVIVOR merge
```
 
 What are some of the options that are available and how many SVs are retained when you run the below code?
```
SURVIVOR merge ~/survivor/inputs/all.txt 250 3 1 1 0 50 ~/survivor/outputs/allVSall.vcf
```
