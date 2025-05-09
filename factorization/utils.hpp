namespace factorization::internal {

template <typename Value, typename Power = Value>
constexpr Value BinPow(Value base, Power power) {
  Value result{1};
  while (power > 0) {
    if (power % 2 != 0) {
      result = result * base;
    }
    base = base * base;
    power /= 2;
  }
  return result;
}

} // namespace factorization::internal
