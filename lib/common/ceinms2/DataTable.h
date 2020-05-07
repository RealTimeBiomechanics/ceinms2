/* -------------------------------------------------------------------------- *
 * CEINMS is a standalone toolbox for neuromusculoskeletal modelling and      *
 * simulation. CEINMS can also be used as a plugin for OpenSim either         *
 * through the OpenSim GUI or API. See https://simtk.org/home/ceinms and the  *
 * NOTICE file for more information. CEINMS development was coordinated       *
 * through Griffith University and supported by the Australian National       *
 * Health and Medical Research Council (NHMRC), the US National Institutes of *
 * Health (NIH), and the European Union Framework Programme 7 (EU FP7). Also  *
 * see the PROJECTS file for more information about the funding projects.     *
 *                                                                            *
 * Copyright (c) 2010-2015 Griffith University and the Contributors           *
 *                                                                            *
 * CEINMS Contributors: C. Pizzolato, M. Reggiani, M. Sartori,                *
 *                      E. Ceseracciu, and D.G. Lloyd                         *
 *                                                                            *
 * Author(s): C. Pizzolato                                                    *
 *                                                                            *
 * CEINMS is licensed under the Apache License, Version 2.0 (the "License").  *
 * You may not use this file except in compliance with the License. You may   *
 * obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.*
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#ifndef ceinms_DataTable_h
#define ceinms_DataTable_h
#include <vector>
#include <string>

#include <ostream>
#include <fstream>


// TODO: insert further checks, operators + and -
namespace ceinms {
template<typename T>
class DataTable;

template<typename T>
std::ostream &operator<<(std::ostream &out, const DataTable<T> &rhs);

template<typename T>
class DataTable {
  public:
    typedef typename std::vector<std::vector<T>>::iterator iterator;
    typedef typename std::vector<std::vector<T>>::const_iterator const_iterator;

    DataTable();
    virtual ~DataTable() {}
    DataTable(size_t nRows, size_t nCols);// do not include time column
    void setLabels(
        const std::vector<std::string> &labels);// do not include time column
    std::vector<std::string> getLabels() const { return labels_; }
    void set(T value, size_t row, size_t col);
    void setTimeColumn(const std::vector<double> &timeColumn);
    void setColumn(size_t col, const std::vector<T> &values);
    void setColumn(const std::string &columnName, const std::vector<T> &values);
    void pushColumn(const std::string &columnName,
        const std::vector<T> &values);
    const T &at(size_t row, size_t col) const;
    T &at(size_t row, size_t col);
    void pushRow(double time, const std::vector<T> &values);
    void setRow(size_t row, const std::vector<T> &values);
    void crop(double startTime, double finalTime);

    static DataTable sum(const DataTable &lhs, const DataTable &rhs);
    static DataTable subtract(const DataTable &lhs, const DataTable &rhs);
    static DataTable multiplyByElement(const DataTable &lhs,
        const DataTable &rhs);
    static DataTable multiplyByScalar(const DataTable &lhs, T scalar);

    std::vector<double> accumulateColumns() const;

    std::string getLabel(size_t col) const { return labels_.at(col); }
    T get(size_t row, size_t col) const { return data_.at(row, col); }
    double getTime(size_t row) const { return time_.at(row); }
    void setTime(double time, std::size_t row);
    T getStartTime() const { return time_.front(); }
    T getFinalTime() const { return time_.back(); }
    std::vector<double> getTimeColumn() const { return time_; }
    const std::vector<T> &getRow(
        size_t row) const;// return const ref is bad.. checking for performance

    std::vector<T> getColumn(size_t col) const;
    std::vector<T> getColumn(const std::string &columnName) const;
    int getColumnIndex(const std::string &columnName) const;
    size_t getNColumns() const { return nCols_; }
    size_t getNRows() const { return nRows_; }

    bool equals(const DataTable &rhs) const;
    friend std::ostream &operator<<<>(std::ostream &output,
        const DataTable &ths);
    void print(const std::string &filename);
    iterator begin() { return data_.begin(); }
    const_iterator begin() const { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator end() const { return data_.end(); }

  private:
    size_t nRows_, nCols_;
    // internal implementation may be changed to something better.. e.g. nice
    // matrices
    std::vector<std::vector<T>> data_;
    std::vector<double> time_;
    std::vector<std::string> labels_;
};
}// namespace ceinms

std::vector<std::string> tokenize(const std::string &line, char sep = ',') {
    std::vector<std::string> result;
    size_t first(0), last(0);
    if (line.find(sep) == std::string::npos)
        result.emplace_back(line);
    else {
        while (last < line.size()) {
            if (line[last] == '\"') {
                last++;
                last = line.substr(last, line.size() - last).find('\"') + last + 1;
            } else if (line[last] == sep) {
                result.emplace_back(line.substr(first, last - first));
                first = ++last;
            } else
                ++last;
        }
        if (last != first) result.emplace_back(line.substr(first, last - first));
    }
    return result;
}

template<typename T>
static ceinms::DataTable<T>
    fillFromCSV(const std::string filename, const char sep = ',', unsigned skiprows = 0) {
    ceinms::DataTable<T> data;
    std::ifstream inF(filename);
    if (!inF.is_open()) {
        std::cout << "Cannot open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    string line;
    for (unsigned i = 0; i < skiprows; ++i)
        getline(inF, line, '\n');
    getline(inF, line, '\n');
    auto labels(tokenize(line, sep));
    data.setLabels(vector<string>(labels.begin() + 1, labels.end()));

    while (!inF.eof()) {
        getline(inF, line);
        if (!line.empty()) {
            vector<string> result(tokenize(line, sep));
            vector<T> values;
            std::transform(std::begin(result),
                std::end(result),
                std::back_inserter(values),
                [](std::string v) { return std::stod(v); });
            data.pushRow(values.at(0), std::vector<T>(values.begin() + 1, values.end()));
        }
    }
    return data;
}

#include "DataTable.cpp"
#endif
