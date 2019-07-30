# Weighted Edge Editing



**Main Algorithm**
```
Input: remaining editing cost k

if k < m_lower_bound(): return no_solution

problem = get_problem(k)
if problem is solved: output edits

if k == 0: return no_solution

solution_found = false
for uv in problem:
	edits := edits u {uv}
	graph.toggle(uv)
	if edit(k - cost(uv)):
		solution_found := true
	graph.toggle(uv)

for uv in problem:
	if uv in edits:
		edits := edits - {uv}

return solution_found

```
