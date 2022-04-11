
# Steps for 1_coordAdjust

## 1. Arrange data into subject directories before ICA analysis

./pregiftdataarrange.sh

## 2. Perform GIFT-ICA analysis manually in Matlab

Toolbox (version 4.0b) can be downloaded here:
https://trendscenter.org/software/gift/

## 3. Adjust original region centers to better match group specific data

Use the coordAdjustGIFT.py script

###Usage

coordAdjustGIFT.py <GIFT-ICA dir> <GIFT-ICA prefix> <NWfile> <Output dir>

Example:
python3 coordAdjustGIFT.py /path/to/diabetes/GIFT-ICA/ diabetes_network rsn_dict.py /path/to/diabetes/GIFT-ICA/
