#include <filesystem>
#include <iostream>
#include <thread>

#include "io.hpp"
#include "mvdst.hpp"

void solve(const std::filesystem::path& inFile,
           const std::filesystem::path& outFile) {
  std::cout << "Now solveing: " << inFile << std::endl;
  graph::Mvdst task{graph::io::read(inFile)};
  //  if (!task.check()) {
  //    std::cout << "HEHMDA  check failed !!!!!!!!!!!!!!!!!!!" << std::endl;
  //  }
  task.approximate();
  graph::io::write(task, outFile);
}

// Чтобы скомпилить, нужно положить библиотеку boost в ../lib/boost/
// python3 imgs.py рисует картинки по файам задач/решений
// Так было выяснено, что распределение вершин совсем не равномерное
int main() {
  //  solve("../Taxicab_64.txt", "../Kurbatov_64.txt");
  //  solve("../Taxicab_128.txt", "../Kurbatov_128.txt");
  std::thread t1(solve, "../Taxicab_512.txt", "../Kurbatov_512.txt");
  //  std::thread t2(solve, "../Taxicab_2048.txt", "../Kurbatov_2048.txt");
  //  solve("../Taxicab_4096.txt", "../Kurbatov_4096.txt");
  t1.join();
  //  t2.join();
  return 0;
}
