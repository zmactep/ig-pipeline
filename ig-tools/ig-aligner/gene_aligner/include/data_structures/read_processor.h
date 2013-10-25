#ifndef READ_PROCESSOR_H_
#define READ_PROCESSOR_H_

#include <omp.h>
#include <iostream>
#include <fstream>
#include "tbb/concurrent_queue.h"
#include "generic_exception.h"

class ReadProcessor {
public:
	ReadProcessor(int num_threads): read(0), num_threads(num_threads){};

	template<typename Reader, typename Handler>
	void readAndProcessSingleThread(Reader& reader, Handler& handler) {
		if(!reader.is_open()) {
			return;
		}

		while (!reader.eof()) {
			Read r;
			reader >> r;
			handler(r);

			if (++read % 1000 == 0) {
			  std::clog << read << " reads processed\r";
			}
		}
	}

	template<typename Reader, typename Handler>
	void readAndProcess(Reader& reader, Handler handler) {
		tbb::concurrent_queue<Read> queue;

		if(!reader.is_open()) {
			return;
		}

		bool stop = false;
#pragma omp parallel shared(reader, queue) firstprivate(handler) num_threads(num_threads)
		{
#pragma omp master
			{
				while (!reader.eof()) {
					Read r;
					reader >> r;
					queue.push(r);
				}
//#pragma omp atomic
				stop = true;
			}

			while(!stop || !queue.empty()) {
				Read r;
				if (queue.try_pop(r)) {
					handler(r);
#pragma omp atomic
					++read;
					if (read % 1000 == 0) {
						std::clog << read << " reads processed\r";
					}
				}
			}
		}
	}

	int get_num_processed_reads() const {
		return read;
	}

private:
	int read;
	int num_threads;
};

#endif /* READ_PROCESSOR_H_ */
