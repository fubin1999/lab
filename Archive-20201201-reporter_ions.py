#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import os
import sys
import json
import re
import time
from progressbar import progressbar # 需安装
from pyteomics import mgf           # 需安装
from collections import namedtuple

#***** 修改下面的参数 *****
REPORTER  = (126.127726, 
			 127.124761, 
			 127.131081, 
			 128.128116, 
			 128.134436, 
			 129.131471, 
			 129.137790, 
			 130.134825, 
			 130.141145, 
			 131.138180)
TOLERANCE = 0.005
#************************

# 常量
SCANS   = 'scans'
CHARGE  = 'charge'
PEPMASS = 'pepmass'

# 储存data的数据结构
Data = namedtuple('Data', 'mz z scan')

# 将工作路径改为当前文件夹
os.chdir(os.path.dirname(sys.argv[0]))


class CSVHandler:
	'''
	处理csv文件的类，负责读取csv文件，提供scan数以供Searcher搜索，
	并将搜索到的结果写入csv文件。
	'''
	def read_data(self, filename):
		data = []
		# 读取csv文件中的m/z、z和scan
		with open(filename, newline='') as csv_file:
			reader = csv.DictReader(csv_file)
			for row in reader:
				mz = row['Observed\nm/z']
				z = row['z']
				# row['Scan #'] 的格式：
				# “controllerType=0 controllerNumber=1 scan=1154”
				scan = row['Scan #'].split('=')[-1]
				data.append(Data(mz, z, scan))
		return data

	def write_data(self, original_file, results, target_path):
		os.makedirs(target_path, exist_ok=True)
		filename = os.path.basename(original_file)[0:-4] + '-result.csv'
		result_path = os.path.join(target_path, filename)
		with open(result_path, 'w', newline='') as fp:
			writer = csv.writer(fp)
			writer.writerows(results)


class Searcher:
	'''从mgf中搜索相应强度的类'''
	def __init__(self, hcd_mgf_file, etd_mgf_file):
		self.hcd_mgf_file = hcd_mgf_file
		self.etd_mgf_file = etd_mgf_file

	def find_reporter(self, spectrum, reporters=REPORTER, tolerance=TOLERANCE):
		'''
		从给定的质谱数据中找到报告离子的强度。
		rpt: 报告离子的mz
		'''
		intens_list = spectrum['intensity array']
		mz_list = spectrum['m/z array']

		# results 是一个列表，列表中的每一项分别对应每一个找到的报告离子的强度
		results = []

		for rpt in reporters:
			mz_max = rpt + tolerance
			mz_min = rpt - tolerance

			# 符合 MIN < mz < MAX 的mz值在列表中的索引组成的列表
			index_list = [mz_list.index(i) for i in mz_list if mz_min < i < mz_max]
			result = max([intens_list[i] for i in index_list]) if index_list else None
			results.append(result)

		return results


	def search(self, scan):
		'''
		根据给定的scan，先在mgf文件里定位母离子数据，
		再查找报告离子的强度。
		:param scan:
		:return result: 报告离子的峰强度
		'''


def get_scan(title):
	'''
	Get scan from the title string.
	title example: "
		20201019SZY-A1.3.3.5 File:"20201019SZY-A1.raw", \
		NativeID:"controllerType=0 controllerNumber=1 scan=3"
		"
	'''
	p = re.compile(r'scan=(\d+)')
	m = p.search(title)
	return m.group(1)


class HCDSearcher(Searcher):

	def __init__(self, hcd_mgf_file, etd_mgf_file):
		super().__init__(hcd_mgf_file, etd_mgf_file)
		self.db = self.generate_db()

	def generate_db(self):
		'''
		db is a dict of dicts. For each item in db, 
		key is the scan number (str),
		value is {
			'm/z array': list,
			'intensity array': list
		}
		'''
		raw_db = mgf.MGF(self.hcd_mgf_file, convert_arrays=0, 
			read_charges=False)

		db = {}
		for spec in raw_db:
			scan = get_scan(spec['params']['title'])
			db[scan] = {
				'm/z array': spec['m/z array'], 
				'intensity array': spec['intensity array']
			}
		return db

	def search(self, scan):
		'''用hcd文件中每个离子的scans，在mfg文件中搜索对应的质谱数据'''
		spectrum = self.db.get(str(scan))
		if spectrum:
			return self.find_reporter(spectrum)
		return [None] * len(REPORTER)


Key = namedtuple('Key', 'pepmass charge')


class ETDSearcher(Searcher):

	def __init__(self, hcd_mgf_file, etd_mgf_file):
		super().__init__(hcd_mgf_file, etd_mgf_file)
		self.hcd_db = self.generate_hcd_db()
		self.etd_db = self.generate_etd_db()

	def generate_hcd_db(self):
		'''
		hcd_db is a dict. For each item in hcd_db:
		* key is a tuple Key(pepmass, charge), in which pepmass is a float, 
		  while charge is an int.
		* value is {
			'm/z array': list,
			'intensity array': list, 
			'scan': str
		}
		'''
		raw_hcd_db = mgf.MGF(self.hcd_mgf_file, convert_arrays=0, 
			read_charges=False)
		hcd_db = {}

		for spec in raw_hcd_db:
			pepmass = spec['params']['pepmass'][0]
			charge = spec['params']['charge'][0].real

			mz_array = spec['m/z array']
			inten_array = spec['intensity array']
			scan = get_scan(spec['params']['title'])

			hcd_db[Key(pepmass, charge)] = {
				'm/z array': mz_array, 
				'intensity array': inten_array, 
				'scan': scan
			}
		return hcd_db

	def generate_etd_db(self):
		'''
		etd_db is a dict. For each item in etd_db:
		* key is the scan number(str).
		* value is a tuple Key(pepmass, charge) (the key in hcd_db)
		'''
		raw_etd_db = mgf.MGF(self.etd_mgf_file, convert_arrays=0, 
			read_charges=False)
		etd_db = {}

		for spec in raw_etd_db:
			scan = get_scan(spec['params']['title'])
			pepmass = spec['params']['pepmass'][0]
			charge = spec['params']['charge'][0].real
			etd_db[scan] = Key(pepmass, charge)
		return etd_db

	def search(self, scan):
		'''用etd文件的每个离子的电荷质荷比，在mgf文件中搜索对应的质谱数据。
		spectrum = self.locate((PEPMASS, CHARGE), (mz, charge))'''
		# 从etd.mgf里通过scan找对应母离子的pepmass和charge
		key = self.etd_db.get(scan)
		# 再从hcd.mgf里通过pepmass和charge找spectrum数据
		spectrum = self.hcd_db.get(key)

		scans_are_similar = abs(int(scan)-int(spectrum['scan'])) < 10
		if spectrum and scans_are_similar:
			return self.find_reporter(spectrum)
		return [None] * len(REPORTER)


class Main:
	'''利用CSVHandler获得的数据，交给Searcher搜索，然后利用CSVHandler将结果写入文件。'''
	def __init__(self, hcd_csv_file, etd_csv_file, hcd_mgf_file, etd_mgf_file, 
			target_path='root'):
		self.hcd_csv_file = hcd_csv_file
		self.etd_csv_file = etd_csv_file
		self.hcd_mgf_file = hcd_mgf_file
		self.etd_mgf_file = etd_mgf_file

		if target_path == 'root':
			root = os.path.dirname(sys.argv[0])
			self.target_path = os.path.join(root, 'result')
		else:
			self.target_path = target_path

	def run(self):
		# 从两个csv文件中读取数据
		handler = CSVHandler()
		self.hcd_data = handler.read_data(self.hcd_csv_file)
		self.etd_data = handler.read_data(self.etd_csv_file)

		# 搜索hcd.csv文件
		hcd_results = []
		hcd_searcher = HCDSearcher(self.hcd_mgf_file, self.etd_mgf_file)
		for hcd_data in self.hcd_data:
			hcd_result = hcd_searcher.search(hcd_data.scan)
			hcd_results.append(hcd_result)

		# 搜索etd.csv文件
		etd_results = []
		etd_searcher = ETDSearcher(self.hcd_mgf_file, self.etd_mgf_file)
		for etd_data in self.etd_data:
			etd_result = etd_searcher.search(etd_data.scan)
			etd_results.append(etd_result)
		
		# 将结果合并至原csv文件中
		handler.write_data(self.hcd_csv_file, hcd_results, self.target_path)
		handler.write_data(self.etd_csv_file, etd_results, self.target_path)



if __name__ == '__main__':
	root = os.path.dirname(sys.argv[0])
	dirs = os.listdir(root)
	dirs = [p for p in dirs if os.path.isdir(os.path.join(root, p))
							if os.path.basename(p) != 'result']

	target_path = os.path.join(root, 'result')

	for folder in progressbar(dirs):
		path = os.path.join(root, folder)
		hcdcsv = os.path.join(path, folder+'-HCD.csv')
		etdcsv = os.path.join(path, folder+'-ETD.csv')
		hcdmgf = os.path.join(path, folder+'-HCD.mgf')
		etdmgf = os.path.join(path, folder+'-ETD.mgf')
		main = Main(hcdcsv, etdcsv, hcdmgf, etdmgf, target_path)
		main.run()

