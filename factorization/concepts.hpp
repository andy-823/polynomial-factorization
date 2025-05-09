#pragma once

#include <concepts>

namespace factorization::internal {

template <typename Field>
concept GaloisField = requires (Field field) {
  { field.Add(Field::Value, Field::Value) } -> std::same_as<typename Field::Value>;
  { field.Negative(Field::Value) } -> std::same_as<typename Field::Value>;
  { field.Multiply(Field::Value, Field::Value) } -> std::same_as<typename Field::Value>;
  { field.Inverse(Field::Value) } -> std::same_as<typename Field::Value>;
};

}  // namespace factorization::internal