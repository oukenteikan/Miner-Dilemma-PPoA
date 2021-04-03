import time
import logging

def init_logger():
    Date = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    BASIC_FORMAT = '%(asctime)s-%(filename)s#%(lineno)d:%(message)s'
    DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(BASIC_FORMAT, DATE_FORMAT)
    chlr = logging.StreamHandler()
    chlr.setFormatter(formatter)
    chlr.setLevel(logging.INFO)
    fhlr = logging.FileHandler(f'./log/log_{Date}.txt')
    fhlr.setFormatter(formatter)
    fhlr.setLevel(logging.INFO)
    logger.addHandler(chlr)
    logger.addHandler(fhlr)
    return logger

if __name__ == "__main__":
    pass