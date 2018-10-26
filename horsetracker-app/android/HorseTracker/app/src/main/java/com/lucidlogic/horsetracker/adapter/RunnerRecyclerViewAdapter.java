package com.lucidlogic.horsetracker.adapter;

import android.databinding.DataBindingUtil;
import android.databinding.ViewDataBinding;
import android.support.v7.widget.RecyclerView;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.TextView;

import com.lucidlogic.horsetracker.BR;
import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.adapter.viewholder.GenericDataBindingViewHolder;
import com.lucidlogic.horsetracker.handler.RaceClickHandler;
import com.lucidlogic.horsetracker.handler.RunnerEntityClickHandler;
import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Runner;

import java.util.ArrayList;
import java.util.List;

import timber.log.Timber;


public class RunnerRecyclerViewAdapter extends RecyclerView.Adapter<GenericDataBindingViewHolder> {

    public List<Runner> getRunners() {
        return runners;
    }

    public void setRunners(List<Runner> runners) {
        this.runners = runners;
    }

    private List<Runner> runners = new ArrayList<>();

    public RunnerRecyclerViewAdapter() {
    }

    public RunnerRecyclerViewAdapter(List<Runner> runners) {
        this.runners = runners;
    }

    @Override
    public GenericDataBindingViewHolder onCreateViewHolder(ViewGroup parent, int viewType) {
        ViewDataBinding viewDataBinding = DataBindingUtil
                .inflate(LayoutInflater.from(parent.getContext()), R.layout.fragment_race, parent, false);
        return new GenericDataBindingViewHolder(viewDataBinding);
    }

    @Override
    public void onBindViewHolder(GenericDataBindingViewHolder holder, int position) {
        Runner postBinding = runners.get(position);
        ViewDataBinding viewDataBinding = holder.getViewDataBinding();
        viewDataBinding.setVariable(BR.runner, postBinding);
        viewDataBinding.setVariable(BR.eventHandler, new RunnerEntityClickHandler());
    }

//    @Override
//    public void onBindViewHolder(final ViewHolder holder, int position) {
//        holder.mItem = runners.get(position);
//        holder.mIdView.setText(runners.get(position).getHorse().getName());
//        holder.mContentView.setText(runners.get(position).getHorse().getProfileUrl());
//
//        holder.mView.setOnClickListener(new View.OnClickListener() {
//            @Override
//            public void onClick(View v) {
//                Timber.i("Nav to horse details");
//            }
//        });
//    }

    @Override
    public int getItemCount() {
        return runners.size();
    }

//    public class ViewHolder extends RecyclerView.ViewHolder {
//        public final View mView;
//        public final TextView mIdView;
//        public final TextView mContentView;
//        public Runner mItem;
//
//        public ViewHolder(View view) {
//            super(view);
//            mView = view;
//            mIdView = (TextView) view.findViewById(R.id.id);
//            mContentView = (TextView) view.findViewById(R.id.content);
//        }
//
//        @Override
//        public String toString() {
//            return super.toString() + " '" + mContentView.getText() + "'";
//        }
//    }
}
