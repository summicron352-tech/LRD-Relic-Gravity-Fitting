from astropy.table import Table
import pandas as pd

def convert_fits_to_csv(fits_filename, csv_filename):
    print(f"正在读取 {fits_filename} ...")
    
    # 1. 使用 astropy 读取 FITS 表格
    # FITS 文件可能会有多个扩展 (extensions)，通常数据在第 1 个扩展中
    data_table = Table.read(fits_filename, format='fits')
    
    # 2. 将数据转换为 Pandas DataFrame（这是一种非常适合导出到 Excel 的格式）
    df = data_table.to_pandas()
    
    # 3. 很多天文 FITS 表格中包含字节类型 (byte) 的字符串，Excel 不太喜欢
    # 我们将其全部解码为普通字符串
    for col in df.select_dtypes([object]).columns:
        df[col] = df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
        
    # 4. 导出为 CSV 文件
    df.to_csv(csv_filename, index=False)
    print(f"✅ 转换成功！数据已保存为: {csv_filename}")
    print(f"总共提取了 {len(df)} 个小红点 (LRDs)。您可以直接用 Excel 打开它了。")

# 执行转换
input_file = 'lrd_table_v1.1.fits'  # 请确保文件名与您下载的一致
output_file = 'Kokorev_LRDs_Full.csv'
convert_fits_to_csv(input_file, output_file)