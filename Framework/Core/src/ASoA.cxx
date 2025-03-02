// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/ASoA.h"
#include "ArrowDebugHelpers.h"
#include "Framework/RuntimeError.h"
#include <arrow/util/key_value_metadata.h>
#include <arrow/util/config.h>

namespace o2::soa
{
void accessingInvalidIndexFor(const char* getter)
{
  throw o2::framework::runtime_error_f("Accessing invalid index for %s", getter);
}
void dereferenceWithWrongType()
{
  throw o2::framework::runtime_error_f("Trying to dereference index with a wrong type in _as<>. Note that if you have several compatible index targets in your process() signature, the last one will be the one actually bound to the getter.");
}
void missingFilterDeclaration(int hash, int ai)
{
  throw o2::framework::runtime_error_f("Null selection for %d (arg %d), missing Filter declaration?", hash, ai);
}

void getterNotFound(const char* targetColumnLabel)
{
  throw o2::framework::runtime_error_f("Getter for \"%s\" not found", targetColumnLabel);
}

void emptyColumnLabel()
{
  throw framework::runtime_error("columnLabel: must not be empty");
}

SelectionVector selectionToVector(gandiva::Selection const& sel)
{
  SelectionVector rows;
  rows.resize(sel->GetNumSlots());
  for (auto i = 0; i < sel->GetNumSlots(); ++i) {
    rows[i] = sel->GetIndex(i);
  }
  return rows;
}

std::shared_ptr<arrow::Table> ArrowHelpers::joinTables(std::vector<std::shared_ptr<arrow::Table>>&& tables)
{
  if (tables.size() == 1) {
    return tables[0];
  }
  for (auto i = 0U; i < tables.size() - 1; ++i) {
    if (tables[i]->num_rows() != tables[i + 1]->num_rows()) {
      throw o2::framework::runtime_error_f("Tables %s and %s have different sizes (%d vs %d) and cannot be joined!",
                                           tables[i]->schema()->metadata()->Get("label").ValueOrDie().c_str(),
                                           tables[i + 1]->schema()->metadata()->Get("label").ValueOrDie().c_str(),
                                           tables[i]->num_rows(),
                                           tables[i + 1]->num_rows());
    }
  }
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns;

  for (auto& t : tables) {
    auto tf = t->fields();
    std::copy(tf.begin(), tf.end(), std::back_inserter(fields));
  }

  auto schema = std::make_shared<arrow::Schema>(fields);

  if (tables[0]->num_rows() != 0) {
    for (auto& t : tables) {
      auto tc = t->columns();
      std::copy(tc.begin(), tc.end(), std::back_inserter(columns));
    }
  }
  return arrow::Table::Make(schema, columns);
}

std::shared_ptr<arrow::Table> ArrowHelpers::concatTables(std::vector<std::shared_ptr<arrow::Table>>&& tables)
{
  if (tables.size() == 1) {
    return tables[0];
  }
  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns;
  assert(tables.size() > 1);
  std::vector<std::shared_ptr<arrow::Field>> resultFields = tables[0]->schema()->fields();
  auto compareFields = [](std::shared_ptr<arrow::Field> const& f1, std::shared_ptr<arrow::Field> const& f2) {
    // Let's do this with stable sorting.
    return (!f1->Equals(f2)) && (f1->name() < f2->name());
  };
  for (size_t i = 1; i < tables.size(); ++i) {
    auto& fields = tables[i]->schema()->fields();
    std::vector<std::shared_ptr<arrow::Field>> intersection;

    std::set_intersection(resultFields.begin(), resultFields.end(),
                          fields.begin(), fields.end(),
                          std::back_inserter(intersection), compareFields);
    resultFields.swap(intersection);
  }

  for (auto& field : resultFields) {
    arrow::ArrayVector chunks;
    for (auto& table : tables) {
      auto ci = table->schema()->GetFieldIndex(field->name());
      if (ci == -1) {
        throw std::runtime_error("Unable to find field " + field->name());
      }
      auto column = table->column(ci);
      auto otherChunks = column->chunks();
      chunks.insert(chunks.end(), otherChunks.begin(), otherChunks.end());
    }
    columns.push_back(std::make_shared<arrow::ChunkedArray>(chunks));
  }

  auto result = arrow::Table::Make(std::make_shared<arrow::Schema>(resultFields), columns);
  return result;
}

arrow::ChunkedArray* getIndexFromLabel(arrow::Table* table, std::string_view label)
{
  auto field = std::find_if(table->schema()->fields().begin(), table->schema()->fields().end(), [&](std::shared_ptr<arrow::Field> const& f) {
    auto caseInsensitiveCompare = [](const std::string_view& str1, const std::string& str2) {
      return std::ranges::equal(
        str1, str2,
        [](char c1, char c2) {
          return std::tolower(static_cast<unsigned char>(c1)) ==
                 std::tolower(static_cast<unsigned char>(c2));
        });
    };

    return caseInsensitiveCompare(label, f->name());
  });
  if (field == table->schema()->fields().end()) {
    o2::framework::throw_error(o2::framework::runtime_error_f("Unable to find column with label %s", label));
  }
  auto index = std::distance(table->schema()->fields().begin(), field);
  return table->column(index).get();
}

void notBoundTable(const char* tableName)
{
  throw o2::framework::runtime_error_f("Index pointing to %s is not bound! Did you subscribe to the table?", tableName);
}

void notFoundColumn(const char* label, const char* key)
{
  throw o2::framework::runtime_error_f(R"(Preslice not valid: table "%s" (or join based on it) does not have column "%s")", label, key);
}

void missingOptionalPreslice(const char* label, const char* key)
{
  throw o2::framework::runtime_error_f(R"(Optional Preslice with missing binding used: table "%s" (or join based on it) does not have column "%s")", label, key);
}

} // namespace o2::soa

namespace o2::framework
{
std::string cutString(std::string&& str)
{
  auto pos = str.find('_');
  if (pos != std::string::npos) {
    str.erase(pos);
  }
  return str;
}

std::string strToUpper(std::string&& str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
  return str;
}
} // namespace o2::framework
