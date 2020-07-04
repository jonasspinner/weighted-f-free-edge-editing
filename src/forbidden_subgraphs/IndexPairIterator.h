
#ifndef WEIGHTED_F_FREE_EDGE_EDITING_INDEXPAIRITERATOR_H
#define WEIGHTED_F_FREE_EDGE_EDITING_INDEXPAIRITERATOR_H


class IndexedVertexPairRange {
public:
    using IndexPair = std::pair<unsigned int, unsigned int>;

private:
    class Iterator {
    private:
        const Vertex *m_vertices{nullptr};
        const IndexPair *m_idx{nullptr};

    public:
        using value_type = VertexPair;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = VertexPair;
        using iterator_category = std::random_access_iterator_tag;

        constexpr Iterator() noexcept {};

        constexpr explicit Iterator(const Vertex *vertices, const IndexPair *pairs) noexcept:
                m_vertices(vertices), m_idx(pairs) {};

        [[nodiscard]] constexpr value_type operator*() const noexcept {
            const auto[i, j] = *m_idx;
            return {m_vertices[i], m_vertices[j]};
        }

        [[nodiscard]] constexpr value_type operator[](std::size_t n) const noexcept {
            return *(*this + n);
        }

        constexpr Iterator &operator++() noexcept {
            ++m_idx;
            return *this;
        }

        [[nodiscard]] constexpr Iterator operator++(int) noexcept {
            Iterator copy(*this);
            ++m_idx;
            return copy;
        }

        constexpr Iterator &operator--() noexcept {
            --m_idx;
            return *this;
        }

        [[nodiscard]] constexpr Iterator operator--(int) noexcept {
            Iterator copy(*this);
            --m_idx;
            return copy;
        }

        constexpr Iterator &operator+=(int n) noexcept {
            m_idx += n;
            return *this;
        }

        constexpr Iterator &operator-=(int n) noexcept {
            m_idx -= n;
            return *this;
        }

        [[nodiscard]] constexpr Iterator operator+(int n) const noexcept {
            Iterator copy(*this);
            copy += n;
            return copy;
        }

        [[nodiscard]] friend constexpr auto operator+(int n, const Iterator &it) noexcept {
            return it + n;
        }

        [[nodiscard]] constexpr Iterator operator-(int n) const noexcept {
            Iterator copy(*this);
            copy -= n;
            return copy;
        }

        [[nodiscard]] constexpr difference_type operator-(const Iterator &other) const noexcept {
            return m_idx - other.m_idx;
        }

        [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
            return m_idx == other.m_idx;
        }

        [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept {
            return !(*this == other);
        }

        [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
            return m_idx < other.m_idx;
        }

        [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
            return m_idx <= other.m_idx;
        }

        [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
            return m_idx > other.m_idx;
        }

        [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
            return m_idx >= other.m_idx;
        }
    };

    const Vertex *m_vertices_begin;
    const IndexPair *m_pairs_begin;
    const IndexPair *m_pairs_end;
public:

    template<class V, class P>
    constexpr explicit IndexedVertexPairRange(V &&vertices, P &&pairs) noexcept:
            m_vertices_begin(std::begin(vertices)), m_pairs_begin(std::begin(pairs)), m_pairs_end(std::end(pairs)) {}

    [[nodiscard]] constexpr auto begin() const noexcept {
        return Iterator{m_vertices_begin, m_pairs_begin};
    }

    [[nodiscard]] constexpr auto end() const noexcept {
        return Iterator{m_vertices_begin, m_pairs_end};
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_INDEXPAIRITERATOR_H
