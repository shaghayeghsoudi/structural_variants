# How to fix smoove errors

In this post I will show how to solve some Smoove issues I came across. I’m using smoove v0.2.8 installed with conda.

## 1. Duphold error

```
atal.nim(49)            sysFatal
Error: unhandled exception: index -1 not in 0 .. 14130 [IndexDefect]

```

*** solution *** The duphold version in the conda smoove environment is v0.2.1, need to update it to version 0.2.3, which is not available on conda.
First download duphold from github and check the installation.



