import random

def sample_discrete_gaussian(std_dev: float, N: int):
    return [SparseEncoder.random_rounding(random.gauss(0, std_dev)) for _ in range(N)]

def sample_zo(p):
    return random.choices([-1, 0, 1], weights=[p/2, 1-p, p/2])[0]

def sample_hwt(h, N: int):
    if N < h:
        raise ValueError("Ring dimension cannot be smaller than the hamming weight")
    empty = N * [0]
    index_to_fill = random.sample(range(N), k=h)
    for index in index_to_fill:
        empty[index] = random.choice([-1, 1])
    return empty

def sample_uniform(ran: int, N: int):
    return [random.randint(-ran/2, ran/2) for _ in range(N)]
