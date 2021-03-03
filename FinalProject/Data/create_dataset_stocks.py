import yfinance 
import pandas as pd
from pandas_datareader import data as pdr

dftickers = pd.read_excel ('NASDAG_CURRENT_INDEX.xlsx')
tickers_list = dftickers['Ticker'].tolist()

data = pdr.get_data_yahoo(tickers_list, start="1985-01-01")
price = data.loc[:, 'Close']
price.to_excel("stockprices_nasdaq.xlsx")