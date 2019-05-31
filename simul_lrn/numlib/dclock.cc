
dclock::dclock(int start_pos, int s_rel, int ende_pos, int e_rel, int order)
{
  reset(start_pos, s_rel, ende_pos, e_rel, order);
}

void dclock::reset(int start_pos, int s_rel, int ende_pos, int e_rel, int order)
{
  start= start_pos;
  start_rel= s_rel;
  ende= ende_pos;
  ende_rel= e_rel;
  ord= order;
  cnt.resize(ord);
  first_time= 1;
  if (ord > 0)
  {
    cnt.set(0, start);
    for (int i= 1; i < ord; i++)
    {
      cnt.set(i, start_rel*cnt[i-1]);
    }
  }
}

int dclock::advance()
{
  bst_set<int> i_set;
  bst_iterator<int> iter(&i_set);
  int forbidden, ibuf;

  if (first_time)
  {
    first_time= 0;
    return 1;
  }
  int level= ord-1;
  while (level > -1)
  {
    cnt.set(level, cnt[level]+1);
    if (((level > 0) && (cnt[level] <= ende_rel*cnt[level-1] + ende))
	|| ((level == 0) && (cnt[level] <= ende)))
    {
      if (level == ord-1)
      {
	i_set.clear();
	for(int j= 0; j < ord; j++)
	{
	  i_set.add(cnt[j]);
	}
	iter.init();
	forbidden= 0;
	while((iter.step()) && (!forbidden))
	{
	  ibuf= iter.c_value();
	  if ((ibuf*(ibuf+1) != 0) && (i_set.in(-(ibuf+1))))
	  {
	    forbidden= 1;
	  }
	}
	if (!forbidden)
	{
	  return 1;
	}
      }
      else
      {
	level++;
	cnt.set(level, start_rel*cnt[level-1]-1);
      }
    }
    else
    {
      level--;
    }
  }
  return 0;
}


int dclock::operator[](int i)
{
  return cnt[i];
}

tnvector<int> dclock::all()
{
  return cnt;
}

tnvector<int> dclock::all_up_ordered()
{
  tnvector<int> vec(cnt);
  int buf;

  for (int i= 0; i < ord; i++)
  {
    for (int j= i+1; j < ord; j++)
    {
      if (vec[i] > vec[j])
      {
	buf= vec[i];
	vec.set(i, vec[j]);
	vec.set(j, buf);
      }
    }
  }
  return vec;
}

tnvector<int> dclock::all_down_ordered()
{
  tnvector<int> vec(cnt);
  int buf;

  for (int i= 0; i < ord; i++)
  {
    for (int j= i+1; j < ord; j++)
    {
      if (vec[i] < vec[j])
      {
	buf= vec[i];
	vec.set(i, vec[j]);
	vec.set(j, buf);
      }
    }
  }
  return vec;
}

