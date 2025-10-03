"""
sudoku_solver.py

Implement the function `solve_sudoku(grid: List[List[int]]) -> List[List[int]]` using a SAT solver from PySAT.
"""

from pysat.formula import CNF
from pysat.solvers import Solver
from typing import List

# TODO: implement encoding and solving using PySAT

def encoder(row: int, col: int, num: int):
    return row * 100 + col * 10 + num

def solve_sudoku(grid: List[List[int]]):
    clauses = CNF()

    # num comes atleast once in each column
    for col in range(1, 10):
        for num in range(1, 10):
            clause = []
            for row in range(1, 10):
                clause.append(encoder(row, col, num))
            clauses.append(clause)
    # for each pair of rows ,num cannot come simultaneously in both rows of same column
            for r1 in range(1, 10):
                for r2 in range(r1 + 1, 10):
                    clauses.append([-encoder(r1, col, num), -encoder(r2, col, num)])

    # there is atmost one number in each cell 
    for row in range(1, 10):
        for col in range(1, 10):
            clause = []
            for num in range(1, 10):
                clause.append(encoder(row, col, num))
            clauses.append(clause)
        #  no two numbers in the same cell
            for n1 in range(1, 10):
                for n2 in range(n1 + 1, 10):
                    clauses.append([-encoder(row, col, n1), -encoder(row, col, n2)])

    # num comes atleast once in each row
        for num in range(1, 10):
            clause = []
            for col in range(1, 10):
                clause.append(encoder(row, col, num))
            clauses.append(clause)

        # num cannot come simultaneously in two columns of same row
            for c1 in range(1, 10):
                for c2 in range(c1 + 1, 10):
                    clauses.append([-encoder(row, c1, num), -encoder(row, c2, num)])

    # Subgrid constraints - for each 3x3 block
    for brow in (1, 4, 7):           
        for bcol in (1, 4, 7):       
            for number in range(1, 10):        
            
                block_positions = []
                for dx in range(3):    
                    for dy in range(3):
                        row = brow + dx
                        col = bcol + dy
                        block_positions.append(encoder(row, col, number))
            
            # At least one cell in the block must contain this number
                clauses.append(block_positions)
            
            # At most one cell in the block can contain this number
                for i in range(len(block_positions)):
                    for j in range(i + 1, len(block_positions)):
                        clauses.append([-block_positions[i], -block_positions[j]])
    # Pre-filled constraint
        for row in range(1, 10):
            for col in range(1, 10):
                value = grid[row - 1][col - 1]
                if value != 0:
                    clauses.append([encoder(row, col, value)])

    # Sudoku solver using PySAT
    with Solver(name='glucose3') as solver:
        for clause in clauses.clauses:     
            solver.add_clause(clause)

        ok = solver.solve()
        model = solver.get_model()

    # Decode solution
    solved = []
    for r in range(9):
        solved.append([0] * 9)

    for literal in model:
        if literal > 0:
            v = literal
            row = v // 100
            col = (v // 10) % 10
            num = v % 10
            solved[row - 1][col - 1] = num

    return solved


