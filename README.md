# Docker_16S
これは16Sの解析用のdockerイメージです。

---

#Docker_16S_pfastqdump
これは16Sの解析用のdockerイメージです / This is a docker image for 16S analysis. (pfastqdump)

# Install & Run

```docker pull petadimensionlab/docker_16s_pfastqdump```

```docker run -it -v /<Users_directory>/:/condir petadimensionlab/docker_16s_pfastqdump```

# Usage
解析対象のデータを<Users_directory>に保存。
Save data to be analyzed in <Users_directory>.

```cd vsearch_ubuntu/```
 
```python3 exec_make_configcsv.py <XXXXX> <yes or no>```

```python3 exec_DLSSRs.py <XXXXX>```

```python3 seq_exec_pfastqdump.py <XXXXX> <yes or no>```

※ <XXXXX> = Folder name of analyzed data
※ <yes or no> = PAIRED is yes , SINGLE is no

Continue analysis to [petadimensionlab/docker_16s_vsearch](https://hub.docker.com/r/petadimensionlab/docker_16s_vsearch/)

---

#
これは16Sの解析用のdockerイメージです / This is a docker image for 16S analysis. (vsearch)

Continued from [petadimensionlab/docker_16s_pfastqdump](https://hub.docker.com/r/petadimensionlab/docker_16s_pfastqdump/).
# Install & Run

```docker pull petadimensionlab/docker_16s_vsearch```

```docker run -it -v /<Users_directory>/:/condir petadimensionlab/docker_16s_vsearch```

# Usage
解析対象のデータを<Users_directory>に保存。
Save data to be analyzed in <Users_directory>.

```cd vsearch_ubuntu/```

```python3 exec_vsearch_full.py <XXXXX> <yes or no>```

※ <XXXXX> = Folder name of analyzed data
※ <yes or no> = PAIRED is yes , SINGLE is no

