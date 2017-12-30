from mpi4py import MPI


class MpiProcess(object):
    def __init__(self, worker, master, data_collector):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.workers_count = self.comm.Get_size() - 1
        self.all_workers = list(range(1, self.workers_count + 1))
        self.free_workers = self.all_workers.copy()
        self.busy_workers = []
        self.worker = worker
        self.master = master
        self.data_collector = data_collector

    def run(self):
        if self.rank == 0:
            self.master(self)
        else:
            self.worker(self)

    def send_to_worker(self, data):
        if not self.free_workers:
            self.collect_result()

        worker_id = self.free_workers.pop()
        self.busy_workers.insert(0, worker_id)

        self.comm.isend(data, dest=worker_id, tag=1)

    def collect_result(self):
        worker_id = self.busy_workers.pop()
        self.free_workers.insert(0, worker_id)

        self._recv_data_from_worker(worker_id)

    def recv_from_master(self):
        data = self.comm.recv(source=0, tag=1)

        return data

    def send_to_master(self, data):
        self.comm.isend(data, dest=0, tag=2)

    def wait_all_workers(self):
        for worker_id in self.busy_workers:
            self._recv_data_from_worker(worker_id)

        for worker_id in self.all_workers:
            self.comm.isend(None, dest=worker_id, tag=1)

    def _recv_data_from_worker(self, worker_id):
        data = self.comm.recv(source=worker_id, tag=2)

        self.data_collector(self, data)
