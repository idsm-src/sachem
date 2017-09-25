#define HEAP_ROW_SIZE           (1024*1024)
#define HEAP_ITEM(heap, idx)    ((heap)->data[(idx) / HEAP_ROW_SIZE][(idx) % HEAP_ROW_SIZE])
#define HEAP_PARENT_IDX(idx)    (((idx) - 1) / 2)
#define HEAP_CHILD_IDX(idx)     ((idx) * 2 + 1)


typedef struct
{
    int32_t id;
    float4 score;
} HeapItem;


typedef struct
{
    HeapItem **data;
    uint32_t size;
    MemoryContext context;
} Heap;


static inline void heap_init(Heap *const heap, int capacity)
{
    int rowCount = ((capacity - 1) / HEAP_ROW_SIZE) + 1;

    heap->data = palloc0(rowCount * sizeof(HeapItem *));
    heap->size = 0;
    heap->context = CurrentMemoryContext;
}


static inline void heap_swap_item(Heap *const heap, int idx1, int idx2)
{
    HeapItem tmp = HEAP_ITEM(heap, idx1);
    HEAP_ITEM(heap, idx1) = HEAP_ITEM(heap, idx2);
    HEAP_ITEM(heap, idx2) = tmp;
}


static inline uint32_t heap_size(Heap *const heap)
{
    return heap->size;
}


static inline HeapItem heap_head(Heap *const heap)
{
    return heap->data[0][0];
}


static inline void heap_add(Heap *const heap, HeapItem const item)
{
    uint32_t idx = heap->size++;
    uint32_t parent = HEAP_PARENT_IDX(idx);

    if(unlikely(heap->data[idx / HEAP_ROW_SIZE] == NULL))
    {
        PG_MEMCONTEXT_BEGIN(heap->context);
        heap->data[idx / HEAP_ROW_SIZE] = palloc(HEAP_ROW_SIZE * sizeof(HeapItem));
        PG_MEMCONTEXT_END();
    }

    HEAP_ITEM(heap, idx) = item;

    while(idx > 0 && item.score > HEAP_ITEM(heap, parent).score)
    {
        heap_swap_item(heap, idx, parent);

        idx = parent;
        parent = HEAP_PARENT_IDX(idx);
    }
}


static inline void heap_remove(Heap *const heap)
{
    heap->size--;

    if(heap->size == 0)
        return;

    HEAP_ITEM(heap, 0) = HEAP_ITEM(heap, heap->size);

    uint32_t idx = 0;
    uint32_t child = HEAP_CHILD_IDX(idx);
    float4 score = HEAP_ITEM(heap, 0).score;

    while(child < heap->size)
    {
        float4 child_score = HEAP_ITEM(heap, child).score;

        if(child + 1 < heap->size)
        {
            float4 right_score = HEAP_ITEM(heap, child + 1).score;

            if(right_score > child_score)
            {
                child_score = right_score;
                child++;
            }
        }

        if(score >= child_score)
            break;

        heap_swap_item(heap, idx, child);

        idx = child;
        child = HEAP_CHILD_IDX(idx);
    }
}
