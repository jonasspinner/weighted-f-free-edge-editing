#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYMATRIX_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYMATRIX_H


namespace dynamic_bitset {
    template<typename Block = unsigned long long>
    class BitRow;

    template<typename Block = unsigned long long>
    class BitMatrix {
        static_assert(std::is_unsigned_v<Block>);
    public:
        using block_type = Block;
        using size_type = Vertex;
        using block_width_type = typename std::vector<Block>::size_type;

        static constexpr block_width_type bits_per_block = std::numeric_limits<Block>::digits;
        static constexpr size_type npos = static_cast<size_type>(-1);

    private:
        size_type m_num_vertices{};
        std::vector<Block> m_blocks{};

        friend
        class BitRow<Block>;

    public:

        class reference {
        private:
            block_type &m_block;
            const block_type m_mask;

            friend class BitMatrix<Block>;

            constexpr reference(block_type &b, block_width_type bit_pos) :
                    m_block(b), m_mask(block_type(1) << bit_pos) {
                assert(bit_pos < bits_per_block);
            }

        public:
            // NOLINTNEXTLINE
            void operator&() = delete;

            [[nodiscard]] constexpr operator bool() const noexcept { // NOLINT(google-explicit-constructor)
                return (m_block & m_mask) != 0;
            }

            [[nodiscard]] constexpr bool operator~() const noexcept {
                return (m_block & m_mask) == 0;
            }

            constexpr reference &flip() noexcept {
                m_block ^= m_mask;
                return *this;
            }

            constexpr reference &operator=(bool x) noexcept {
                if (x) {
                    m_block |= m_mask;
                } else {
                    m_block &= ~m_mask;
                }
                return *this;
            }

            constexpr reference &operator=(const reference &rhs) noexcept {
                *this = bool(rhs);
                return *this;
            }

            constexpr reference &operator|=(bool x) noexcept {
                if (x) {
                    m_block |= m_mask;
                }
                return *this;
            }

            constexpr reference &operator&=(bool x) noexcept {
                if (!x) {
                    m_block &= ~m_mask;
                }
                return *this;
            }

            constexpr reference &operator^=(bool x) noexcept {
                if (x) {
                    m_block ^= m_mask;
                }
                return *this;
            }

            constexpr reference &operator-=(bool x) noexcept {
                if (x) {
                    m_block &= ~m_mask;
                }
                return *this;
            }
        };

        using const_reference = bool;

        constexpr BitMatrix() noexcept = default;

        explicit constexpr BitMatrix(size_type num_vertices) :
                m_num_vertices(num_vertices), m_blocks(num_vertices * calc_num_blocks_per_row(m_num_vertices)) {}

        BitMatrix(const BitMatrix &b) = default;

        BitMatrix &operator=(const BitMatrix &b) = default;

        BitMatrix(BitMatrix &&src) noexcept: m_num_vertices(src.m_num_vertices),
                                             m_blocks(std::move(src.m_blocks)) {
            src.m_num_vertices = 0;
        };

        BitMatrix &operator=(BitMatrix &&src) noexcept {
            m_num_vertices = src.m_num_vertices;
            m_blocks = std::move(src.m_blocks);
            src.m_num_vertices = 0;
        };

        ~BitMatrix() {
            assert(check_invariants());
        };

        void swap(BitMatrix &other) noexcept {
            using std::swap;
            swap(m_blocks, other.m_blocks);
            swap(m_num_vertices, other.m_num_vertices);
        }

        constexpr void clear() noexcept {
            m_blocks.clear();
            m_num_vertices = 0;
        }

        class ConstRowView {
            const block_type *const m_blocks;
            const size_type m_num_vertices{};

            friend class BitMatrix<Block>;

            friend class BitRow<Block>;

            constexpr ConstRowView(const block_type *blocks, size_type num_vertices) :
                    m_blocks(blocks), m_num_vertices(num_vertices) {}

        public:
            ConstRowView() = delete;

            ConstRowView(const ConstRowView &other) = delete;

            // NOLINTNEXTLINE
            void operator&() = delete;

            ConstRowView &operator=(const ConstRowView &rhs) = delete;

            [[nodiscard]] constexpr bool test(size_type pos) const {
                assert(pos < m_num_vertices);
                return (m_blocks[block_index(pos)] & bit_mask(pos)) != 0;
            }

            [[nodiscard]] constexpr const_reference operator[](size_type pos) const {
                assert(pos < m_num_vertices);
                return test(pos);
            }

            [[nodiscard]] constexpr size_type size() const noexcept {
                return m_num_vertices;
            }

            class ConstBlockView {
                const block_type *const m_blocks;
                const size_type m_num_blocks;
            public:
                ConstBlockView(const block_type *blocks, size_type num_blocks) : m_blocks(blocks),
                                                                                 m_num_blocks(num_blocks) {}

                [[nodiscard]] constexpr auto begin() const noexcept {
                    return m_blocks;
                }

                [[nodiscard]] constexpr auto end() const noexcept {
                    return m_blocks + m_num_blocks;
                }

                [[nodiscard]] constexpr auto size() const noexcept {
                    return m_num_blocks;
                }
            };

            [[nodiscard]] constexpr auto blocks() const noexcept {
                return ConstBlockView{m_blocks, calc_num_blocks_per_row(m_num_vertices)};
            }

            friend std::ostream &operator<<(std::ostream &os, const ConstRowView &view) {
                for (Vertex u = 0; u < view.m_num_vertices; ++u) {
                    os << view[u] << " ";
                }
                return os;
            }

            class IndexRange {
                class IndexIterator {
                private:
                    const block_type *m_blocks{nullptr};
                    size_type m_num_vertices{};
                    size_type m_pos{};

                public:
                    using value_type = size_type;
                    using difference_type = std::ptrdiff_t;
                    using pointer = void;
                    using reference = value_type;
                    using iterator_category = std::forward_iterator_tag;

                    constexpr IndexIterator() noexcept = default;

                    constexpr IndexIterator(const block_type *blocks, size_type num_vertices, size_type pos = 0) :
                            m_blocks(blocks), m_num_vertices(num_vertices), m_pos(pos) {
                        assert(m_pos <= m_num_vertices);
                        if (m_pos < m_num_vertices && !(m_blocks[block_index(pos)] & bit_mask(pos))) {
                            ++*this;
                        }
                    };

                    [[nodiscard]] constexpr value_type operator*() const noexcept {
                        assert(m_pos < m_num_vertices);
                        return m_pos;
                    }

                    IndexIterator &operator++() noexcept {
                        assert(m_pos < m_num_vertices);
                        const auto blocks = m_blocks;
                        const auto num_blocks = calc_num_blocks_per_row(m_num_vertices);
                        const auto num_bits = m_num_vertices;

                        ++m_pos;

                        if (m_pos >= num_bits) {
                            m_pos = num_bits;
                            return *this;
                        }

                        auto block_idx = block_index(m_pos);
                        const auto idx = bit_index(m_pos);

                        assert(block_idx < num_blocks);
                        const auto remaining_block = blocks[block_idx] >> idx;
                        if (remaining_block != Block(0)) {
                            m_pos += static_cast<size_type>(detail::ctz(remaining_block));
                            assert(m_pos <= num_bits);
                            return *this;
                        }

                        ++block_idx;
                        while (block_idx < num_blocks && (blocks[block_idx] == Block(0)))
                            ++block_idx;

                        m_pos = block_idx >= num_blocks
                                ? num_bits
                                : block_idx * bits_per_block + static_cast<size_type>(detail::ctz(blocks[block_idx]));

                        assert(m_pos <= num_bits);
                        return *this;
                    }

                    [[nodiscard]] IndexIterator operator++(int) noexcept { // NOLINT(cert-dcl21-cpp)
                        IndexIterator copy(*this);
                        ++*this;
                        return copy;
                    }

                    [[nodiscard]] constexpr bool operator==(const IndexIterator &other) const noexcept {
                        return m_pos == other.m_pos;
                    }

                    [[nodiscard]] constexpr bool operator!=(const IndexIterator &other) const noexcept {
                        return !(*this == other);
                    }

                    [[nodiscard]] constexpr bool operator<(const IndexIterator &other) const noexcept {
                        return m_pos < other.m_pos;
                    }

                    [[nodiscard]] constexpr bool operator<=(const IndexIterator &other) const noexcept {
                        return m_pos <= other.m_pos;
                    }

                    [[nodiscard]] constexpr bool operator>(const IndexIterator &other) const noexcept {
                        return m_pos > other.m_pos;
                    }

                    [[nodiscard]] constexpr bool operator>=(const IndexIterator &other) const noexcept {
                        return m_pos >= other.m_pos;
                    }
                };

                const block_type *m_blocks;
                size_type m_num_vertices;
            public:
                constexpr IndexRange(const block_type *blocks, size_type num_vertices) :
                        m_blocks(blocks), m_num_vertices(num_vertices) {}

                [[nodiscard]] constexpr auto begin() const {
                    return IndexIterator{m_blocks, m_num_vertices};
                }

                [[nodiscard]] constexpr auto end() const {
                    return IndexIterator{m_blocks, m_num_vertices, m_num_vertices};
                }
            };

            [[nodiscard]] constexpr auto indices() const {
                return IndexRange{m_blocks, m_num_vertices};
            }
        };

        class RowView {
            block_type *const m_blocks;
            const size_type m_num_vertices{};
            const size_type m_num_blocks{};

            friend class BitMatrix<Block>;

            friend class BitRow<Block>;

            RowView(block_type *blocks, size_type num_vertices) : m_blocks(blocks), m_num_vertices(num_vertices),
                                                                  m_num_blocks(calc_num_blocks_per_row(num_vertices)) {}

        public:
            RowView() = delete;

            RowView(const RowView &other) = delete;

            // NOLINTNEXTLINE
            void operator&() = delete;

            constexpr RowView &operator=(const ConstRowView &rhs) {
                assert(size() == rhs.size());
                std::copy(rhs.m_blocks, std::next(rhs.m_blocks, num_blocks()), m_blocks);
                return *this;
            }

            constexpr RowView &operator=(const RowView &rhs) {
                if (&rhs == this)
                    return *this;
                *this = static_cast<ConstRowView>(rhs);
                return *this;
            }

            constexpr void swap(RowView &rhs) noexcept {
                using std::swap;
                assert(size() == rhs.size());
                assert(check_non_overlapping(m_blocks, rhs.m_blocks, num_blocks()));
                std::swap(m_blocks, std::next(m_blocks, num_blocks()), rhs.m_blocks);
            }

            [[nodiscard]] constexpr operator ConstRowView() const {  // NOLINT(google-explicit-constructor)
                return ConstRowView(m_blocks, m_num_vertices);
            }

            constexpr RowView &operator&=(const ConstRowView &rhs) {
                assert(size() == rhs.size());
                assert(check_non_overlapping(m_blocks, rhs.m_blocks, num_blocks()));
                for (size_type i = 0; i < num_blocks(); ++i)
                    m_blocks[i] &= rhs.m_blocks[i];
                return *this;
            }

            constexpr RowView &operator|=(const ConstRowView &rhs) {
                assert(size() == rhs.size());
                assert(check_non_overlapping(m_blocks, rhs.m_blocks, num_blocks()));
                for (size_type i = 0; i < num_blocks(); ++i)
                    m_blocks[i] |= rhs.m_blocks[i];
                return *this;
            }

            constexpr RowView &operator^=(const ConstRowView &rhs) {
                assert(size() == rhs.size());
                assert(check_non_overlapping(m_blocks, rhs.m_blocks, num_blocks()));
                for (size_type i = 0; i < num_blocks(); ++i)
                    m_blocks[i] ^= rhs.m_blocks[i];
                return *this;
            }

            constexpr RowView &operator-=(const ConstRowView &rhs) {
                assert(size() == rhs.size());
                assert(check_non_overlapping(m_blocks, rhs.m_blocks, num_blocks()));
                for (size_type i = 0; i < num_blocks(); ++i)
                    m_blocks[i] &= ~rhs.m_blocks[i];
                return *this;
            }

            constexpr RowView &set() noexcept {
                std::fill(m_blocks, std::next(m_blocks, num_blocks()), ~static_cast<Block>(0));
                zero_unused_bits();
                return *this;
            }

            constexpr RowView &set(size_type pos) {
                assert(pos < m_num_vertices);
                m_blocks[block_index(pos)] |= bit_mask(pos);
                return *this;
            }

            constexpr RowView &reset() noexcept {
                std::fill(m_blocks, std::next(m_blocks, num_blocks()), Block(0));
                return *this;
            }

            constexpr RowView &reset(size_type pos) {
                assert(pos < m_num_vertices);
                m_blocks[block_index(pos)] &= ~bit_mask(pos);
                return *this;
            }

            constexpr RowView &flip() noexcept {
                for (size_type i = 0; i < num_blocks(); ++i)
                    m_blocks[i] = ~m_blocks[i];
                zero_unused_bits();
                return *this;
            }

            constexpr RowView &flip(size_type pos) {
                assert(pos < m_num_vertices);
                m_blocks[block_index(pos)] ^= bit_mask(pos);
                return *this;
            }

            [[nodiscard]] constexpr bool test(size_type pos) const {
                assert(pos < m_num_vertices);
                return (m_blocks[block_index(pos)] & bit_mask(pos)) != 0;
            }

            [[nodiscard]] constexpr reference operator[](size_type pos) {
                assert(pos < m_num_vertices);
                return reference(m_blocks[block_index(pos)], bit_index(pos));
            }

            [[nodiscard]] constexpr const_reference operator[](size_type pos) const {
                assert(pos < m_num_vertices);
                return test(pos);
            }

            [[nodiscard]] constexpr size_type size() const noexcept {
                return m_num_vertices;
            }

            [[nodiscard]] constexpr size_type num_blocks() const noexcept {
                return m_num_blocks;
            }

        private:
            constexpr void zero_unused_bits() noexcept {
                assert(num_blocks() == calc_num_blocks_per_row(m_num_vertices));

                // if != 0 this is the number of bits used in the last block
                const block_width_type extra_bits = count_extra_bits();

                if (extra_bits != 0)
                    highest_block() &= (Block(1) << extra_bits) - 1;
            }

            [[nodiscard]] constexpr block_width_type count_extra_bits() const noexcept {
                return bit_index(size());
            }

            [[nodiscard]] constexpr Block &highest_block() {
                assert(size() > 0 && num_blocks() > 0);
                return m_blocks[num_blocks() - 1];
            }

            [[nodiscard]] constexpr const Block &highest_block() const {
                assert(size() > 0 && num_blocks() > 0);
                return m_blocks[num_blocks() - 1];
            }

            [[nodiscard]] constexpr bool check_invariants() const noexcept {
                const block_width_type extra_bits = count_extra_bits();
                if (extra_bits > 0) {
                    const block_type mask = (~block_type(0)) << extra_bits;
                    if ((highest_block() & mask) != 0)
                        return false;
                }
                if (num_blocks() != calc_num_blocks_per_row(size()))
                    return false;

                return true;
            };

            [[nodiscard]] static constexpr bool
            check_non_overlapping(const Block *first_row, const Block *second_row, size_type num_blocks) noexcept {
                return first_row + num_blocks <= second_row || second_row + num_blocks <= first_row;
            };

        };

        [[nodiscard]] constexpr RowView mut_row(Vertex u) {
            const auto idx = u * calc_num_blocks_per_row(m_num_vertices);
            return RowView(&m_blocks[idx], m_num_vertices);
        }

        [[nodiscard]] constexpr ConstRowView row(Vertex u) const {
            const auto idx = u * calc_num_blocks_per_row(m_num_vertices);
            return ConstRowView(&m_blocks[idx], m_num_vertices);
        }

        [[nodiscard]] constexpr reference operator[](VertexPair uv) {
            const auto[u, v] = uv;
            assert(u < m_num_vertices);
            assert(v < m_num_vertices);
            const auto idx = u * calc_num_blocks_per_row(m_num_vertices) + block_index(u);
            assert(idx < m_blocks.size());
            return reference(m_blocks[idx], bit_index(v));
        }

        [[nodiscard]] constexpr bool test(VertexPair uv) const {
            const auto[u, v] = uv;
            assert(u < m_num_vertices);
            assert(v < m_num_vertices);
            const auto idx = u * calc_num_blocks_per_row(m_num_vertices) + block_index(u);
            assert(idx < m_blocks.size());
            return (m_blocks[idx] & bit_mask(v)) != 0;
        }

        [[nodiscard]] constexpr bool operator[](VertexPair uv) const {
            return test(uv);
        }

        friend std::ostream &operator<<(std::ostream &os, const BitMatrix &matrix) {
            for (Vertex u = 0; u < matrix.m_num_vertices; ++u) {
                os << matrix.row(u) << "\n";
            }
            return os;
        }

    private:
        [[nodiscard]] constexpr bool check_invariants() noexcept {
            for (Vertex u = 0; u < m_num_vertices; ++u) {
                if (!this->mut_row(u).check_invariants()) {
                    return false;
                }
            }
            return true;
        };

        [[nodiscard]] constexpr static size_type block_index(size_type pos) noexcept {
            return pos / bits_per_block;
        }

        [[nodiscard]] constexpr static block_width_type bit_index(size_type pos) noexcept {
            return static_cast<block_width_type>(pos % bits_per_block);
        }

        [[nodiscard]] constexpr static Block bit_mask(size_type pos) noexcept {
            return Block(1) << bit_index(pos);
        }

        [[nodiscard]] constexpr static size_type calc_num_blocks_per_row(size_type num_bits) noexcept {
            return num_bits / bits_per_block
                   + static_cast<size_type>(num_bits % bits_per_block != 0);
        }
    };

    template<typename Block>
    inline void swap(BitMatrix<Block> &lhs, BitMatrix<Block> &rhs) noexcept(noexcept(lhs.swap(rhs))) {
        lhs.swap(rhs);
    }

    template<class Block>
    class BitRow {
        using size_type = typename BitMatrix<Block>::size_type;
        std::vector<Block> m_blocks{};
        size_type m_num_vertices;

        using RowView = typename BitMatrix<Block>::RowView;
        using ConstRowView = typename BitMatrix<Block>::ConstRowView;
        using reference = typename BitMatrix<Block>::reference;
        using const_reference = typename BitMatrix<Block>::const_reference;
    public:
        BitRow() = delete;

        explicit constexpr BitRow(size_type num_vertices) :
                m_blocks(BitMatrix<Block>::calc_num_blocks_per_row(num_vertices)), m_num_vertices(num_vertices) {}

        constexpr BitRow &operator=(const BitRow &rhs) = default;

        constexpr BitRow(BitRow &&other) noexcept:
                m_blocks(std::move(other.m_blocks)), m_num_vertices(other.m_num_vertices) {
            other.m_num_vertices = 0;
        }

        constexpr BitRow &operator=(BitRow &&rhs) noexcept {
            m_num_vertices = rhs.m_num_vertices;
            m_blocks = std::move(rhs.m_blocks);
            rhs.m_num_vertices = 0;
            return *this;
        }

        void swap(BitRow &other) {
            using std::swap;
            swap(m_blocks, other.m_blocks);
            swap(m_num_vertices, other.m_num_vertices);
        }

    private:
        [[nodiscard]] constexpr RowView mut_view() {
            return RowView(m_blocks.data(), m_num_vertices);
        }

        [[nodiscard]] constexpr ConstRowView view() const {
            return ConstRowView(m_blocks.data(), m_num_vertices);
        }

    public:

        constexpr BitRow &operator=(const ConstRowView &rhs) {
            mut_view() = rhs;
            return *this;
        }

        constexpr void swap(RowView &rhs) noexcept {
            mut_view().swap(rhs);
        }

        constexpr BitRow &operator&=(const ConstRowView &rhs) {
            mut_view() &= rhs;
            return *this;
        }

        constexpr BitRow &operator|=(const ConstRowView &rhs) {
            mut_view() |= rhs;
            return *this;
        }

        constexpr BitRow &operator^=(const ConstRowView &rhs) {
            mut_view() ^= rhs;
            return *this;
        }

        constexpr BitRow &operator-=(const ConstRowView &rhs) {
            mut_view() -= rhs;
            return *this;
        }

        constexpr BitRow &set() noexcept {
            mut_view().set();
            return *this;
        }

        constexpr BitRow &set(size_type pos) {
            mut_view().set(pos);
            return *this;
        }

        constexpr BitRow &reset() noexcept {
            mut_view().reset();
            return *this;
        }

        constexpr BitRow &reset(size_type pos) {
            mut_view().reset(pos);
            return *this;
        }

        constexpr BitRow &flip() noexcept {
            mut_view().flip();
            return *this;
        }

        constexpr BitRow &flip(size_type pos) {
            mut_view().flip();
            return *this;
        }

        [[nodiscard]] constexpr bool test(size_type pos) const {
            return view().test(pos);
        }

        [[nodiscard]] constexpr reference operator[](size_type pos) {
            return mut_view()[pos];
        }

        [[nodiscard]] constexpr const_reference operator[](size_type pos) const {
            return view()[pos];
        }

        [[nodiscard]] constexpr size_type size() const noexcept {
            return m_num_vertices;
        }

        [[nodiscard]] constexpr size_type num_blocks() const noexcept {
            return BitMatrix<Block>::calc_num_blocks_per_row(m_num_vertices);
        }

        friend std::ostream &operator<<(std::ostream &os, const BitRow &row) {
            return os << row.view();
        }

        [[nodiscard]] constexpr auto indices() const {
            return view().indices();
        }
    };

    template<typename Block>
    inline void swap(BitRow<Block> &lhs, BitRow<Block> &rhs) noexcept(noexcept(lhs.swap(rhs))) {
        lhs.swap(rhs);
    }
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYMATRIX_H
