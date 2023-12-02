from get_testcase_data import get_testcase_data
from create_ics import create_ics
from run_model import run_model
from check_results import check_results

nCells = 10242
n1 = 16
n2 = 32

get_testcase_data()

create_ics(nCells)

run_model(nCells, n1)

run_model(nCells, n2)

check_results(n1, n2)
