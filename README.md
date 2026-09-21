# Outline

# Environment

Python 3.13</br>
encoding UTF-8

# Usage
コマンドプロンプト、ターミナル等で[`main.py`](main.py)を実行する。
```
$ python main.py
```

解析対象のファイル、結果の出力先はmain関数内の引数で指定する。
```
if __name__=='__main__':

    main('example/thrust.csv', 'output/sample')
```
第1引数：解析対象のファイル、第2引数：結果の出力先

# Set up
必要なPythonのライブラリのインストール
```
$ python -m pip install -r requirements.txt
```