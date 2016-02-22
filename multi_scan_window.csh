mkdir $1

date > granite.4.log
python scan.py granite.4a --w 11.5,13.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.4b --w 66.5,68.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.4c --w 121.5,123.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.4d --w 176.5,178.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.4e --w 231.5,233.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.4f --w 286.5,288.5,11.5,13.5 --sxi 0.1 --syi 0.1 --m --s

date > granite.5.log
python scan.py granite.5a --w 11.5,13.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.5b --w 66.5,68.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.5c --w 121.5,123.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.5d --w 176.5,178.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.5e --w 231.5,233.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.5f --w 286.5,288.5,66.5,68.5 --sxi 0.1 --syi 0.1 --m --s

date > granite.6.log
python scan.py granite.6a --w 11.5,13.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.6b --w 66.5,68.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.6c --w 121.5,123.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.6d --w 176.5,178.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.6e --w 231.5,233.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.6f --w 286.5,288.5,121.5,123.5 --sxi 0.1 --syi 0.1 --m --s

date > granite.7.log
python scan.py granite.7a --w 11.5,13.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.7b --w 66.5,68.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.7c --w 121.5,123.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.7d --w 176.5,178.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.7e --w 231.5,233.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.7f --w 286.5,288.5,176.5,178.5 --sxi 0.1 --syi 0.1 --m --s

date > granite.8.log
python scan.py granite.8a --w 11.5,13.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.8b --w 66.5,68.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.8c --w 121.5,123.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.8d --w 176.5,178.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.8e --w 231.5,233.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.8f --w 286.5,288.5,231.5,233.5 --sxi 0.1 --syi 0.1 --m --s

date > granite.9.log
python scan.py granite.9a --w 11.5,13.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.9b --w 66.5,68.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.9c --w 121.5,123.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.9d --w 176.5,178.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.9e --w 231.5,233.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s
python scan.py granite.9f --w 286.5,288.5,286.5,288.5 --sxi 0.1 --syi 0.1 --m --s

mv granite.* $1
